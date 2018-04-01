# -*- coding: utf-8 -*-
"""
Prog: nt_shannon_entropy
Desc: Generate and visualize Shannon-log2 entropy from nucleotide data.

Req. Std. Dependencies: argparse, math, sys
Opt. Std Dependencies: os--for file searching
Req. External Dependencies: Bio.SeqIO, Bio.Alphabet.IUPAC, 
Bio.SeqRecord.SeqRecord, Bio.Seq.Seq, numpy, matplotlib.pyplot

Author: Michael J. Bale (michael.bale@nih.gov)
Version: 1.0.2
Date: 03-31-2018
"""

import argparse
import math
import sys
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
import matplotlib.pyplot as plt

def edit_Seq(seqs):
    """Changes all '-' in sequence set to 'N'"""
    new_seqs = list()
    for record in seqs:
        new_seqs.append(SeqRecord(Seq(str(record.seq).replace('-','N'),
                                      IUPAC.ambiguous_dna),
                                  id=record.id, description = ""))
    #next record
    return new_seqs
#/def edit_Seq

def _get_nucleotide_composition(nt):
    """determine what nucleotides are in supplied nt
       most important in ambiguous nucleotide-containing sequences:
       E.g. V = A, G, and C    
    """
    one_way_nt = ['C', 'A', 'T', 'G', '-']
    two_way_ambiguity_avec = ['Mac', 'Rag', 'Wat', 'Scg', 'Yct', 'Ktg']
    three_way_ambiguity_avec = ['Vagc', 'Hatc', 'Datg', 'Btgc']
    two_way_ambiguity_sans = ['M', 'R', 'W', 'S', 'Y', 'K']
    three_way_ambiguity_sans = ['V', 'H', 'D', 'B']
    if nt.upper() in one_way_nt:
        return [one_way_nt.index(nt.upper())]
    elif nt.upper() in two_way_ambiguity_sans:
        return [one_way_nt.index(two_way_ambiguity_avec[
                two_way_ambiguity_sans.index(nt.upper())][1].upper()),
                one_way_nt.index(two_way_ambiguity_avec[
                two_way_ambiguity_sans.index(nt.upper())][2].upper())]
    elif nt.upper() in three_way_ambiguity_sans:
        return [one_way_nt.index(three_way_ambiguity_avec[
                three_way_ambiguity_sans.index(nt.upper())][1].upper()),
                one_way_nt.index(three_way_ambiguity_avec[
                three_way_ambiguity_sans.index(nt.upper())][2].upper()),
                one_way_nt.index(three_way_ambiguity_avec[
                three_way_ambiguity_sans.index(nt.upper())][3].upper())]
    else:
        return [0,1,2,3]
    #/if
#/def _get_nucleotide_composition

def generate_allele_frequencies(SequenceSet, handle_indels = False):
    """generate frequencies of nucleotide composition at each position in sequence set"""
    allele_frequencies = np.zeros(5*len(SequenceSet[0].seq)).reshape(
            5,len(SequenceSet[0].seq))
    for index in range(0, len(SequenceSet)):
        for nt in range(0, len(SequenceSet[0].seq)):
            nuc_value = _get_nucleotide_composition(
                    str(SequenceSet[index].seq)[nt])
            for val in nuc_value:
                allele_frequencies[val,nt]+=1/len(nuc_value)
            #next val
        #next nt
    #next index
    if not (handle_indels == True):
        allele_frequencies[4,:] = 0
    #/if
    norm_allele_frequencies = allele_frequencies/np.sum(
                allele_frequencies,axis = 0)
    return np.round(norm_allele_frequencies, decimals = 5)
#/def generate_allele_frequencies

def calc_Shannon(FreqData):
    """Calculate Shannon-Log2 Entropy data"""
    shannon_array = np.zeros(FreqData.shape[1])
    #For each column: -Sum(column)*log2(column)
    for column in range(0,FreqData.shape[1]):
        for val in FreqData[:,column]:
            if val == 0:
                next
            else:
                shannon_array[column]-= val * math.log2(val)
            #/if
        #next val
    #next column
    return np.round(np.transpose(shannon_array),decimals = 5)
#/def calc_Shannon

def print_Shannon(Shannon, fasta_file, threshold = 0.6):
    """Print Shannon data in graphical and return single-column text file"""
    #setup path strings for output
    shannon_output_path = fasta_file + "_ShannonEntropy_values.txt"
    shannon_graph_output_path = fasta_file + "_ShannonEntropy_Graph.png"
    
    #save shannon data as single-column text file
    np.savetxt(shannon_output_path, Shannon, fmt = '%.5f', delimiter = '\t',
               header='''
               Output of Shannon entropy of sequence set where the value in 
               column j is given by -sum(p_i * log_2(p_i)) where p_i is the 
               frequency of nucleotide i in sequence position j.
               ''')
    #Create Graph
    fig, ax = plt.subplots(figsize = (11,8.5))
    ax.invert_yaxis()
    greater_than_threshold = [i for i, val in enumerate(Shannon) if val >
                              threshold]    
    ax.scatter(range(0,len(Shannon)), Shannon)
    ax.scatter(greater_than_threshold, Shannon[greater_than_threshold], 
               color = 'r')
    fig.savefig(shannon_graph_output_path)
#/def print_Shannon
    
def main(*argv):
    parser = argparse.ArgumentParser(
        description='''
                        Program to generate and plot Shannon-log2 Entropy of 
                        Nucleotide Data. IUPAC Ambiguous coding is supported 
                        for data, but trailing indels as '-' may not be 
                        used--instead use an N.
                        ''',
        epilog=
        '''
                        Any questions/issues/bugs may be reported to
                        Michael J. Bale at michael.bale@nih.gov
                        ''', add_help=True
                        )
    parser.add_argument('--fasta', '-f', required=True,
                        help='''
                        Input file containing DNA/RNA sequences. This file 
                        should be aligned and not have an outgroup. IUPAC
                        notations are legal.
                        ''', nargs=1, type=str
                       )
    parser.add_argument('--Threshold', '-t', required=False, type=float, 
                        nargs = 1, default=0.6,
                        help='''
                        Threshold for Shannon Plot. Default is set to 0.6. All
                        bit values at or above this value will be visualized as
                        red instead of blue on plot.
                        '''
                        )
    parser.add_argument('--MissingData', '-M', nargs='?', const='N', default='-',
                        help='''
                        If indels are meaningful in data, use this flag to
                        handle them. Default is to handle them as nothing,
                        addition of the flag without paramter handles all '-'
                        as 'N'. 'Real' handles '-' as a spacer for insertions
                        or deletions in the raw sequence data.
                        '''
                       )    
    args = vars(parser.parse_args())
    
    #handle input values and exit on error
    
    try:
        sequence_set = list(SeqIO.parse(args['fasta'][0],"fasta"))
    except IOError:
        sys.exit("Error opening sequence file: " + args['fasta'][0])
    #/try

    if not((args['MissingData'].upper() == 'REAL') and (
            args['MissingData'].upper() == 'N') and (
            args['MissingData'].upper() == '-')):
        sys.exit("Invalid MissingData argument: " + args['MissingData'] + 
                 "\n\nPlease Use 'Real', 'N', or '-'")
    #/if
        
    if args['MissingData'] == 'N':
        sequence_set = edit_Seq(sequence_set)
    #/if
    
    handle_indels = False
    if args['MissingData'].upper() == 'REAL':
        handle_indels = True
    #/if
    
    ###########################################################################
    
    #calculate allele frequencies and calculate and print shannon values
    allele_freq = generate_allele_frequencies(sequence_set, handle_indels)
    shannon_data = calc_Shannon(allele_freq)
    print_Shannon(shannon_data, args['fasta'][0], args['Threshold'])
#/def main
    
    
#if used as script    
if __name__ == "__main__":
    main(sys.argv[1:])