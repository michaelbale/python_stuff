#! /mnt/nasapps/development/python/3.4.9/bin/python3

# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 17:20:10 2020

@author: Michael
"""

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
import math
import random
import re
import sys
import tqdm
from Bio import SeqIO
import multiprocessing as mp
import numpy as np

def hammingdist(s1, s2):
    hd = 0
    
    if len(s1) != len(s2):
        raise ValueError("Undefined for strings of different lengths")
    
    for c1,c2 in zip(s1, s2):
        if c1 != '-' and c2 != '-' and c1 != c2:
            hd += 1
    return math.log(hd + 1)

def getHDMat(sequences):
    mat = [
            [hammingdist(seq1.seq, seq2.seq) for seq2 in sequences]
            for seq1 in sequences
          ]
    return mat
    
def calcK(pops, taxa_distance, group_list):
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
    K = 0
    groups = [tup[0] for tup in pops]
    sizes = np.array([tup[1] for tup in pops])
    w = (sizes-2)/(len(group_list) - 2*len(pops))
    combsizes = sizes*(sizes - 1) / 2
    k = [0]*len(w)
    for i in range(len(group_list)):
        for j in range(i+1, len(group_list)):
            if group_list[i] == group_list[j]:
                k[groups.index(group_list[i])] += taxa_distance[i][j]
                
    K = np.dot(w, k/combsizes)
    return K
#/def calc_correl_coeff
    
def calcK_bootstrap(pop_sizes, taxa_distance, group_list):
    random.shuffle(group_list)
    return calcK(pop_sizes, taxa_distance, group_list)

def panmixLoop(group_file, taxa_distance, taxa, num_reps):
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
    
    pop_sizes = [(x, group_list.count(x)) for x in set(group_list)]
    kst = calcK(pop_sizes, taxa_distance, group_list)
    k_dist = [calcK_bootstrap(pop_sizes, taxa_distance, group_list) for i in tqdm.tqdm(range(num_reps))]
    
#    for i in tqdm.tqdm(range(num_reps)):
#        random.shuffle(group_list)
#        k_dist.append(calcK(pop_sizes, taxa_distance, group_list))
#    #next i
    p_val = sum(1 for i in k_dist if i <= kst)/num_reps

    return (kst, p_val)
#/def analyze_homochronous_gruops


    
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
                        Input File containing phylogenetic tree.
                        ''', nargs=1, type=str
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
    #/end parser

    args = vars(parser.parse_args(argv))

    group_file = args['group_file']
    num_reps = args['num_bootstraps']

    seq_file_name = args['in_file'][0]

    sequences = list(SeqIO.parse(seq_file_name, "fasta"))
    taxa = [seq.id for seq in sequences]
    HD_mat = getHDMat(sequences)


    

    kst, p_val = panmixLoop(group_file, HD_mat, taxa, num_reps)
    
    print("Kst = " + str(round(kst, 3)))
    print("P-value:" + str(p_val))
    print("Num bootstraps used was " + str(num_reps))
#/def main

if __name__ == "__main__":
    main(sys.argv[1:])
#/end if
#/end tree_correl_coeff.py
