# -*- coding: utf-8 -*-
"""
@author: Michael J. Bale
"""
from Bio import AlignIO
from Bio import Phylo
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np
import argparse
import os
import subprocess
import sys
#N.B. Change below as needed
sys.path.append("C:/Users/Michael/Desktop/python_code/")
#sys.path.append("C:/Users/balemj/Desktop/Scripts/Python/")
import root_to_tip_analysis as rtt

def _run_dir_analysis(in_files, num_files, time_input, root, out_prefix, 
                      store_dists):
    p_vals = list()
    f_stat = list()
    slopes = list()
    for file in in_files:
        phylip_aln_convert = rtt.convert_alignment(file)
        path_to_sms = "/home/balemj/RScripts/SMS_PhyML/sms-1.8.1/sms.sh"
        subprocess_args = [path_to_sms, '-i', phylip_aln_convert, '-d', 'nt',
                           '-t', '-s', 'SPR']
        subprocess.run(subprocess_args, stdout=open('root_to_tip.log',"w"))
        tree_file_name = phylip_aln_convert + "_phyml_tree.txt"
        tree = Phylo.read(open(tree_file_name, "rU"))
        tree.root_with_outgroup({'name': root})
        paired_taxa_distance = rtt.get_root_to_tip_distances(tree, root)
        agg_results =  rtt.analyze_longitudinal_groups(time_input,
                                                       paired_taxa_distance)
        os.remove(phylip_aln_convert)
        slopes.append(agg_results[0])
        f_stat.append(agg_results[1])
        p_vals.append(agg_results[2])
        if store_dists:
            out_file = file + "_R2TAnalysis.txt"
            file_handle = open(out_file, "w")
            file_handle.write(rtt._generate_header(agg_results))
            for key, val in agg_results[3].items():
                file_handle.write('{}\t{}\t{}\n'.format(key,
                                  str(np.around(float(val[0]), 3)),
                                  str(np.around(val[1], 8))))
    
    out_file = out_prefix + "_windowResults.txt"
    fh = open(out_file, "w")
    fh.write('FileName\tR2TSlope\tF_statistic\tSlope_pVal\n')
    for i in range(num_files):
        fh.write('{}\t{}\t{}\t{}\n'.format(in_files[i], 
                 str(np.around(float(slopes[i], 3))),
                 str(np.around(float(f_stat[i], 3))),
                 str(np.around(float(p_vals[i], 3))))
                 )
        

def _run_list_analysis(input_arg, alignment, root, out_prefix, time_input,
                       store_dists):
    seq_list = list(SeqIO.parse(alignment, "fasta"))
    p_vals = list()
    f_stat = list()
    slopes = list()
    with open(input_arg, "r") as i_fh:
        lines = (line.strip() for line in i_fh)
        lines = list(line for line in lines if line)
    for window in lines:
        new_aln = list()
        for seq in seq_list:
            new_seq = "".join([seq.seq[i] for i in window])
            new_rec = SeqRecord(new_seq, seq.id, description="")
            new_aln.append(new_rec)
            
        tmp_out = alignment + "_tmp.phy"
        _ = AlignIO.write(new_aln, tmp_out, "phylip-relaxed")
        path_to_sms = "/home/balemj/RScripts/SMS_PhyML/sms-1.8.1/sms.sh"
        subprocess_args = [path_to_sms, '-i', tmp_out, '-d', 'nt',
                           '-t', '-s', 'SPR']
        subprocess.run(subprocess_args, stdout=open('root_to_tip.log',"w"))
        tree_file_name = tmp_out + "_phyml_tree.txt"
        tree = Phylo.read(open(tree_file_name, "rU"))
        tree.root_with_outgroup({'name': root})
        paired_taxa_distance = rtt.get_root_to_tip_distances(tree, root)
        agg_results =  rtt.analyze_longitudinal_groups(time_input,
                                                       paired_taxa_distance)
        os.remove(tmp_out)
        slopes.append(agg_results[0])
        f_stat.append(agg_results[1])
        p_vals.append(agg_results[2])
        if store_dists:
            out_file = input_arg + "_R2TAnalysis.txt"
            file_handle = open(out_file, "w")
            file_handle.write(rtt._generate_header(agg_results))
            for key, val in agg_results[3].items():
                file_handle.write('{}\t{}\t{}\n'.format(key,
                                  str(np.around(float(val[0]), 3)),
                                  str(np.around(val[1], 8))))
    
    out_file = out_prefix + "_windowResults.txt"
    fh = open(out_file, "w")
    fh.write('FileName\tR2TSlope\tF_statistic\tSlope_pVal\n')
    for i in range(lines):
        fh.write('{}\t{}\t{}\t{}\n'.format("Window_1", 
                 str(np.around(float(slopes[i], 3))),
                 str(np.around(float(f_stat[i], 3))),
                 str(np.around(float(p_vals[i], 3))))
                )
                    
def main(argv):
    parser = argparse.ArgumentParser(
        prog='window_analysis.py',
        description='''
                        Analysis of slices of sequence alignments for the
                        detection of selective pressures on subsequences
                        within a global alignment.
                        ''',
        epilog='''
                        Any questions/issues/bugs may be reported to Michael J.
                        Bale at michael.bale@nih.gov
                        ''',
        add_help=True
        )
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--is_dir', '-d', default=False, action='store_true',
                        help='''
                        Determines type of input. If is_dir, takes all fasta
                        files within the supplied input directory from --input
                        and runs the indicated analysis.
                        ''')
    group.add_argument('--is_fasta', '-f', default=True, action='store_false',
                       help='''
                       If singular input file is a fasta file, use this flag
                       and parsing will proceed with a fasta file using for 
                       the analysis.
                       ''')
    group.add_argument('--alignment', '-a', default=None, type=str,
                       help='''
                       Supply fasta file with for the coordinates supplied in
                       --input flag.
                       ''')
    parser.add_argument('--input', '-i', required=True, type=str, nargs=1,
                        help='''
                        Input File or directory for analysis. Supplied
                        argument will be handled as indicated by the --is_dir
                        flag. If --is_dir is present, --input will be handled
                        as a directory, else --is_dir is not present and
                        --input will be handled as a fasta file.
                        ''')
    parser.add_argument('--root', '-r', required=True,
                        help='''
                       Taxon name that is to be the outgroup in the tree to be
                       generated. Used in argument passing to 
                       root_to_tip_analysis.py commands.
                       ''', nargs=1, type=str)
    parser.add_argument('--time_file', '-t', required=True,
                        help='''
                       File with discriminating sequences types and dates or
                       groups as tsv format. Used in argument passing to 
                       root_to_tip_analysis.py commands.
                       ''')
    parser.add_argument('--out_file', '-o',
                        help='''
                        File prefix for output values if desired.
                        ''')
    parser.add_argument('--store_dists', '-s', default=False,
                        action='store_true',
                        help='''
                        Stores root-to-tip distances for each call of run -- 
                        recommended for continuous windows, not necessarily
                        recommended for random windows. Only applicable if
                        --is_fasta is not chosen. --is_fasta run automatically
                        stores distances in separate file. See documentation
                        on root_to_tip_analysis.py.
                        ''')
#    parser.add_argument('--do_panmixia', '-p', default=False,
#                        action='store_true',
#                        help='''
#                        Option to include an estimation of K* subdivision
#                        heterogeneity estimation by the method of Hudson,
#                        Kaplan, and Boos.
#                        ''')
    args = vars(parser.parse_args(argv))
    is_dir = args['is_dir']
    is_fasta = args['is_fasta']
    alignment = args['alignment']
    input_arg = args['input'][0]
    root = args['root'][0]
    time_input = args['time_file']
    out_prefix = args['out_file']
    store_dists = args['store_dists']
    assertion_error_message = "Input directory \"" + input_arg + "\" does not exist"
    if is_dir:
        assert os.path.exists(input_arg), assertion_error_message
        in_files = [file for file in os.listdir(input_arg) if 'fas' in file[-3]]
        num_files = len(in_files)
        _run_dir_analysis(in_files, num_files, time_input, root, out_prefix,
                          store_dists)
    elif (is_fasta) and (alignment is None):
        print("Calling \"root_to_tip_analysis.py\" on supplied arguments")
        rtt.main(['-i', input_arg, '-r', root, '-t', time_input])
    else:
        _run_list_analysis(input_arg, alignment, root, out_prefix, time_input,
                           store_dists)
    
if __name__ == "__main__":
    main(sys.argv[1:])
#/end if
#/end window_analysis.py