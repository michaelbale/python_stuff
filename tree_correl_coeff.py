"""
Program: tree_correl_coeff.py

Description: Pairwise branch-length distances from maximum likelihood tree.
Tree is generated from best-fit model as defined by the Smart Model Selection
by Lefort, Jongeuville, and Gascuel MBE 2017 and estimated in PhyML from
Guidon et al. Sys Biol 2010.

Required Standard Dependencies: random, re, subprocess, sys
Required Non-Standard Dependencies: argparse
Required External Modules: Bio.Phylo, Bio.AlignIO, numpy
Required External Programs: Smart Model Selection, PhyML
(both are freely-downloadable at: www.atgc-montpellier.fr)

@version: 1.0.1
@author: Michael J. Bale (michael.bale@nih.gov)
Version Date: 15 May, 2018
"""

import argparse
import random
import re
import subprocess
import sys
from Bio import AlignIO
from Bio import Phylo
import numpy as np





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

def get_pairwise_distance(tree, taxa):
    """
    Create list of lists where list_ij is the branch length distance between
    taxa[i] and taxa[j] on the provided tree. If i == j, value is 0.0. These
    values are not counted in the final calculations.


    Keyword args:

        tree -- Bio.Phylo tree object
        taxa -- list of all taxa in tree

    return:

        paired_distance -- list of list
    """
    paired_distance = list()
    tmp = list()
    for taxon1 in taxa:
        for taxon2 in taxa:
            if taxon1 == taxon2:
                tmp.append(0.0)
            else:
                tmp.append(tree.distance(taxon1, taxon2))
            #/end if
        #next taxon2
        paired_distance.append(tmp)
        tmp = list()
    #next taxon1
    return paired_distance
#/def get_pairwise_distance

def calc_correl_coeff(taxa, taxa_distance, group_list):
    """
    Calculation of intersequence branch length distance correlation
    coefficient.


    Keyword args:

        taxa -- list of all taxa
        taxa_distance -- list of lists containing pairwise taxa distances
        group_list -- list containing group designations; to be shuffled for
                      bootstrapping.

    Internal args:

        x_list -- numerical designation of intersequence compartment relation.
                  x[i] = 0 if same compartment; else 1
        y_list -- single list of intersequence distances paired with x_list

    Internal calculations:

        r_num = sum(sum(x[i] - x_average) * sum(y[i] - y_average))
        r_denom = sqrt(sum((x[i] - x_average)**2)*sum((y[i]-y_average)**2))
        correlation coefficient = r_num/r_denom

    return:

        correlation coefficient
    """
    x_list = list()
    y_list = list()
    for i in range(0, len(taxa)):
        for j in range(0, len(taxa)):
            if not i == j:
                y_list.append(np.around(taxa_distance[i][j], 8))
                tmp = 0 if group_list[i] == group_list[j] else 1
                x_list.append(tmp)
            #/end if
        #next j
    #next i

    x_bar = np.average(x_list)
    y_bar = np.average(y_list)
    r_num = np.sum((x_list-x_bar)*(y_list-y_bar))
    r_denom = np.sqrt(np.sum((x_list-x_bar)**2) * np.sum((y_list-y_bar)**2))
    return r_num/r_denom
#/def calc_correl_coeff

def analyze_node_distance(group_file, taxa_distance, taxa, num_reps):
    """
    Open supplied time file and pass to @method calc_correl_coeff. Uses
    bootstrapped correlation coefficients to determine significance.


    Keyword args:

        group_file -- file containing group designations for taxa list
        taxa_distance -- paired intersequence branch length distances
        taxa -- list of taxa names
        num_reps -- number of bootstrap replicates

    Internal args:

        r -- standard correlation coefficient from @method calc_correl_coeff
        r_dist -- bootstrapped values of correlation coeffiencts
        r* -- number of r_dist[i] that satisfy r_dist[i]>=r

    Internal calculations:

        p_value = r*/num_reps

    return:

        p_value
    """
    with open(group_file, "r") as in_file:
        file_content = in_file.readlines()
    #/end with
    in_file.close()
    file_content = [x.strip() for x in file_content]
    file_dict = {x.split()[0] : x.split()[1] for x in file_content}
    group_list = [0] * len(taxa)
    for pattern, group in file_dict.items():
        for i in range(0, len(taxa)):
            if re.search(pattern, taxa[i]) is not None:
                group_list[i] = group
            #/end if
        #next key, value
    #next pattern, group
    if not len(group_list) == len(taxa_distance):
        raise Exception('time_file does not cover all sequences')
    #/end if
    correl = calc_correl_coeff(taxa, taxa_distance, group_list)
    r_dist = list()
    for i in range(0, num_reps):
        random.shuffle(group_list)
        r_dist.append(calc_correl_coeff(taxa, taxa_distance, group_list))
    #next i
    p_val = sum(1 for i in r_dist if i >= correl)/num_reps

    return (correl, p_val)
#/def analyze_homochronous_gruops

def _generate_header(correl, p_val):
    """Generates header for output file"""
    header = "Analysis of Intersequence distances by permutation test agianst"
    header += "Correlation Coefficient. Intersequence distances are printed"
    header += " below in tsv-matrix form.\n"
    header += "Correlation is X = 0 if seq1 and seq2 are in same group,"
    header += " X = 1 otherwise. Y = total branch length between seq1 and "
    header += " seq2\n"
    header += "Correlation Coefficient: {:+.3E}".format(correl)
    header += "(p: {:.4E})\n".format(np.around(p_val, 4))
    header += "\nIntersequence Distance Matrix\n\n\t"
    return header
#/def _generate_header

def main(argv):
    """Main method for @Prog tree_correl_coeff.py"""

    #begin parser
    parser = argparse.ArgumentParser(
        prog='tree_correl_coeff.py',
        description='''
                    Analysis of Root To Tip distances from maximum
                    likelihood tree. Tree is generated from best-fit 
                    model as defined by the Smart Model Selection by
                    Lefort, Jongeuville, and Gascuel MBE 2017 and
                    estimated in PhyML from Guidon et al. Sys Biol
                    2010. This program analyzes data from
                    homochronous samples by permutation of the GCC
                    using pairwise branch length distance as described
                    in Critchlow et al. Mathematical and Computer
                    Modeling 2000. 
                    ''',
        epilog='''
                    Any questions/issues/bugs may be reported to Michael J.
                    Bale at michael.bale@nih.gov
                ''',
        add_help=True
        )
    parser.add_argument('--in_file', '-i', required=True,
                        help='''
                        Input File containing DNA/RNA sequences.
                        ''', nargs=1, type=str
                       )
    parser.add_argument('--style', '-s', default="fasta",
                        help='''
                        File type for input file. Default is fasta. Input file
                        will be modified to phylip-relaxed for use in the SMS
                        and PhyML programs.
                        '''
                       )
    parser.add_argument('--group_file', '-g', required=True,
                        help='''
                        File with discriminating sequences types and dates or
                        groups as tsv format.
                        '''
                       )
    parser.add_argument('--num_bootstraps', '-n', default=10000, type=int,
                        help='''
                        Number of bootstrap permutations to run for
                        significance testing of the calculation correlation
                        coefficient.
                        '''
                       )
    parser.add_argument('--out_file', '-o',
                        help='''
                        File name for root to tip values and associated
                        analyses if applicable.
                        '''
                       )
    #/end parser

    args = vars(parser.parse_args(argv))
    if args['out_file'] is None:
        out_file = args['in_file'][0] + "_R2TAnalysis.txt"
    else:
        out_file = args['out_file']
    #/end if

    group_file = args['group_file']
    num_reps = args['num_bootstraps']

    #N.B. Change 'path_to_sms' to correct path on individual machine
    path_to_sms = "/home/balemj/RScripts/SMS_PhyML/sms-1.8.1/sms.sh"
    print("Creating Phylip File...")
    phylip_aln_file = convert_alignment(args['in_file'][0], args['style'])
    print("Phylip File Generated...\nRunning SMS and PhyML...")
    subprocess_args = [path_to_sms, "-i", phylip_aln_file, "-d", "nt", "-t"]
    subprocess.run(subprocess_args, stdout=open("root_to_tip.log", "w"))
    print("PhyML tree created...")
    tree_file_name = phylip_aln_file + "_phyml_tree.txt"

    tree = Phylo.read(tree_file_name, "newick")
    taxa = [clade.name for clade in tree.get_terminals()]
    taxa_distance = get_pairwise_distance(tree, taxa)


    file_handle = open(out_file, "w")
    print("Analyzing Root To Tip Distances...")

    correl, p_val = analyze_node_distance(group_file, taxa_distance,
                                          taxa, num_reps)

    file_handle.write(_generate_header(correl, p_val))
    print("Writing Results...")
    for i in range(0, len(taxa)):
        file_handle.write(taxa[i] + '\t')
    #next i
    file_handle.write('\n')

    for i in range(0, len(taxa)):
        file_handle.write(taxa[i] + '\t')
        for j in range(0, len(taxa)):
            file_handle.write(str(np.around(taxa_distance[i][j], 8)) + '\t')
        #next j
        file_handle.write('\n')
    #next i
    file_handle.close()
    print("Runtime complete")
    print("Delete intermediary Files if desired")
#/def main

if __name__ == "__main__":
    main(sys.argv[1:])
#/end if
#/end tree_correl_coeff.py
