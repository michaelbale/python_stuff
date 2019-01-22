"""
Program: root_to_tip_analysis.py

Description: Analysis of Root To Tip distances from maximum likelihood tree.
Tree is generated from best-fit model as defined by the Smart Model Selection
by Lefort, Jongeuville, and Gascuel MBE 2017 and estimated in PhyML from
Guidon et al. Sys Biol 2010. This program only analyzes data from longitudinal
samples.

Required Standard Dependencies: re, subprocess, sys
Required Non-Standard Dependencies: argparse
Required External Modules: Bio.Phylo, Bio.AlignIO, numpy, scipy.stats
Required External Programs: Smart Model Selection, PhyML
(both are freely-downloadable at: www.atgc-montpellier.fr)

@version: 1.0.1
@author: Michael J. Bale (michael.bale@nih.gov)
Version Date: 10 May, 2018
"""

import argparse
import re
import subprocess
import sys
from Bio import AlignIO
from Bio import Phylo
import numpy as np
import scipy.stats

def convert_alignment(in_file, style="fasta"):
    """
    Converts input alignment file to relaxed phylip format for SMS/PhyML.


    Keyword args:

        in_file -- file name for input alignment file to be converted
        style="fasta" -- file type of input file--default is fasta

    return:

        out_file -- file name of output phylip-relaxed file
    """
    out_file = in_file + "_R2Taln.phy"
    msa = AlignIO.read(open(in_file), style)
    #N.B. @var __ is unused output of AlignIO.write() function
    __ = AlignIO.write(msa, open(out_file, "w"), "phylip-relaxed")
    return out_file
#/def convert_alignment

def get_root_to_tip_distances(tree, root):
    """
    Creates dict of taxa : root to tip distance pairs.


    Keyword args:

        tree -- Bio.Phylo rooted tree object
        root -- string containing name of root of @param tree

    return:

        {taxa:root_to_tip_distance}
    """
    taxa = [clade.name for clade in tree.get_terminals()]
    del taxa[taxa.index(root)]
    root_to_tip_distances = [tree.distance(taxon, root) for taxon in taxa]
    return dict(zip(taxa, root_to_tip_distances))
#/def get_root_to_tip_distances

def linear_regression(group_dist_dict):
    """
    Performs Ordinary Least Squares regression on supplied data. Significance
    testing is performed by F-test against inferred linear model.


    Keyword args:

        group_dist_dict -- dictionary with entries taxon : (time-point, R2T)

    Internal args:

        x values -- time-point from group_dist_dict value tuple
        y values -- R2T from group_dist_dict value tuple

    Internal Calculations:

        slope -- (sum(x*y)-1/len(x)*sum(x)*sum(y)) /
                 (sum(x**2)-1/len(x)*sum(x)**2)
        intercept -- y_average - slope * x_average
        model residuals -- sum((y_model - y_average)**2)
        squared error = sum((y_model - y)**2)
        model DoF = 1
        DoF = len(x) - 2
        F statistic = model_residuals*DoF/(model_DoF*squared_error)
        p value = 1-F_distribution_cumulative_density(F, model_DoF, DoF)

    return:

        list(slope, f_statistic, p_value, group_dist_dict)
    """
    x_list = list()
    y_list = list()
    for xy_pair in group_dist_dict.values():
        x_list.append(float(xy_pair[0]))
        y_list.append(float(xy_pair[1]))
    #next xy_pair
    x_array = np.array(x_list)
    y_array = np.array(y_list)

    x_transform = np.vstack([x_array, np.ones(len(x_array))]).T
    model = np.linalg.lstsq(x_transform, y_array)
    slope, intercept = model[0]
    y_bar = np.average(y_array)
    model_residuals = np.sum((np.array([slope*x_array + intercept])-y_bar)**2)
    squared_error = model[1][0]
    f_statistic = (model_residuals*(len(x_list)-2))/squared_error
    p_value = 1 - scipy.stats.f.cdf(f_statistic, 1, len(x_list)-2)

    return list([slope, f_statistic, p_value, group_dist_dict])

def analyze_longitudinal_groups(time_file, paired_taxa_distance):
    """
    Open supplied time file to group R2T dists and pass to
    @method linear_regression().
    """
    with open(time_file, "r") as in_file:
        file_content = in_file.readlines()
    #/end with
    in_file.close()
    file_content = [x.strip() for x in file_content]
    file_dict = {x.split()[0] : x.split()[1] for x in file_content}
    try:
        all(isinstance(float(v), float) for v in file_dict.values())
    except ValueError as err:
        print(str(err))
        quit()
    #/end try-except
    group_dist_dict = dict()
    for pattern, group in file_dict.items():
        for key, value in paired_taxa_distance.items():
            if re.search(pattern, key) is not None:
                group_dist_dict[key] = tuple([group, value])
            #/end if
        #next key, value
    #next pattern, group
    if not len(group_dist_dict) == len(paired_taxa_distance):
        raise Exception('time_file does not cover all sequences')
    #/end if
    results = linear_regression(group_dist_dict)
    return results
#/def analyze_longitudinal_groups

def _generate_header(long_results):
    """Generates header for output file"""
    header = "Analysis of Root to Tip distances by F-test agianst"
    header += "Linear Regression.\nX-Y pairs are printed below in"
    header += "form SeqName\\tTime point\\tRoot To Tip Distance.\n"
    header += "Slope of Regression: {:.3E} ".format(long_results[0])
    header += "(F: {}; ".format(str(np.around(long_results[1], 3)))
    header += "p: {:.3E})\n".format(long_results[2])
    header += "\nSeqName\tTime Point\tRoot To Tip Distance\n"
    return header
#/def _generate_header

def main(argv):
    """Main method for @Prog root_to_tip_analysis.py"""
    parser = argparse.ArgumentParser(
        prog='lon_root_to_tip_analysis.py',
        description='''
                        Analysis of Root To Tip distances from maximum
                        likelihood tree. Tree is generated from best-fit
                        model as defined by the Smart Model Selection by
                        Lefort, Jongeuville, and Gascuel MBE 2017 and
                        estimated in PhyML from Guidon et al. Sys Biol
                        2010. This program only analyzes data from
                        longitudinal samples.
                        ''',
        epilog='''
                        Any questions/issues/bugs may be reported to Michael J.
                        Bale at michael.bale@nih.gov
                        ''',
        add_help=True
        )
    parser.add_argument('--in_file', '-i', required=True,
                        help='''
                        Input File containing DNA/RNA sequences. File should
                        have an outgroup (see Option '--root').
                        ''', nargs=1, type=str
                       )
    parser.add_argument('--root', '-r', required=True,
                        help='''
                       Taxon name that is to be the outgroup in the tree to be
                       generated.
                       ''', nargs=1, type=str)
    parser.add_argument('--style', '-s', default="fasta",
                        help='''
                        File type for input file. Default is fasta. Input file
                        will be modified to phylip-relaxed for use in the SMS
                        and PhyML programs.
                        ''')
    parser.add_argument('--time_file', '-t', required=True,
                        help='''
                       File with discriminating sequences types and dates or
                       groups as tsv format.
                       ''')
    parser.add_argument('--out_file', '-o',
                        help='''
                        File name for root to tip values and associated
                        analyses if applicable.
                        ''')
    args = vars(parser.parse_args(argv))
    if args['out_file'] is None:
        out_file = args['in_file'][0] + "_R2TAnalysis.txt"
    else:
        out_file = args['out_file']
    #/end if
    path_to_sms = "/home/balemj/RScripts/SMS_PhyML/sms-1.8.1/sms.sh"
    root = args['root'][0]
    print("Creating Phylip File...")
    phylip_aln_file = convert_alignment(args['in_file'][0], args['style'])
    print("Running SMS and PhyML...")
    subprocess_args = path_to_sms + " -i " + phylip_aln_file + " -d nt -t -s SPR"
    subprocess.run(subprocess_args, stdout=open("root_to_tip.log", "w"))
    print("PhyML tree created...")
    tree_file_name = phylip_aln_file + "_phyml_tree.txt"

    tree = Phylo.read(open(tree_file_name, "rU"), "newick")
    tree.root_with_outgroup({'name': root})
    paired_taxa_distance = get_root_to_tip_distances(tree, root)
    file_handle = open(out_file, "w")
    print("Analyzing Root To Tip Distances...")
    long_results = analyze_longitudinal_groups(args['time_file'],
                                               paired_taxa_distance)
    print("Writing Results...")
    file_handle.write(_generate_header(long_results))
    for key, val in long_results[3].items():
        file_handle.write('{}\t{}\t{}\n'.format(key,
                                                str(np.around(float(val[0]),
                                                              3)),
                                                str(np.around(val[1], 8))))
    #/next k,v
    file_handle.close()
    print("Runtime complete")
    print("Delete intermediary Files if desired")
#/def main

if __name__ == "__main__":
    main(sys.argv[1:])
#/end if
#/end root_to_tip_analysis.py
