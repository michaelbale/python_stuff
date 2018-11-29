# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 14:39:57 2018

@author: Michael
"""
import argparse
import csv
import itertools
import math
import os
import random as rand
import re
import sys
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def _random_gen(low, high):
    while True:
        yield rand.randrange(low, high)

def _store_as_fasta(alignment, windows, random, len_windows,
                    num_windows, step, in_file):
    if random:
        new_dir_name = in_file + "_randWindows_dir"
    else:
        new_dir_name = in_file + "_contWindows_dir"
    file_name = re.split(r"[\\/]", in_file)[-1]
    if not os.path.isdir(new_dir_name):
        os.makedirs(new_dir_name)
    new_in_file = new_dir_name + "\\" + file_name
    if(random):
        _store_as_array_text(windows, random, len_windows,
                             num_windows, step, new_in_file)
    for i in range(0, num_windows):
        window = windows[i]
        new_aln = list()
        for sequence in alignment:
            seq = str(sequence.seq)
            tmp_list = list()
            tmp_list.extend([seq[index] for index in window])
            new_seq = "".join(tmp_list)
            new_aln.append(SeqRecord(Seq(new_seq, IUPAC.IUPACUnambiguousDNA), sequence.id, description=""))
        new_file_name = new_dir_name + "\\" + file_name + "-window"
        new_file_name += str(i) + "_Size" + str(len_windows)
        if not random:
            new_file_name += "_Step" + str(step)
        new_file_name += ".fas"
        SeqIO.write(new_aln, new_file_name, "fasta")
    
def _store_as_array_text(windows, random, len_windows,
                         num_windows, step, in_file):
    if random:
        file_name = in_file + "_randWindows-size" + str(len_windows) + "_"
        file_name += str(num_windows) + "windows_asArrays.txt"
    else:
        file_name = in_file + "_contWindows-size" + str(len_windows) + "_"
        file_name += str(num_windows) + "windows_step"
        file_name += str(step) + "_windows_asArrays.txt"
    with open(file_name, 'w') as fh:
        wr = csv.writer(fh, quoting=csv.QUOTE_MINIMAL)
        for wind in windows:
            wr.writerow(wind)
    fh.close()


def generate_random_windows(length, len_window, N):
    windows = list()
    for _ in range(N):
        generator = _random_gen(0, length)
        window = set()
        for x in itertools.takewhile(lambda x: len(window) < len_window,
                                     generator):
            window.add(x)
        window = list(window)
        windows.append(window)
    return windows

def generate_contig_windows(len_window, N, step):
    windows = list()
    for wind in range(N):
        start_pos = wind * step
        end_pos = start_pos + len_window
        window = [x for x in range(start_pos, end_pos)]
        windows.append(window)
    return windows

def main(argv):

    parser = argparse.ArgumentParser(
        prog='make_windows_from_fasta.py',
        description='''
                            Creates slices of sequences alignments as 
                            contiguous or random windows. Can store these
                            windows as individual fasta files or as a file
                            contatining the coordinates for each window.
                            ''',
        epilog='''
                            Any questions/issues/bugs/ may be reported to
                            Michael J. Bale at michael.bale@nih.gov
                            ''',
        add_help=True
        )
    parser.add_argument('--in_file', '-i', required=True,
                        help='''
                        Input File containing sequences to generate windows
                        from.
                        ''', nargs=1, type=str
                       )
    parser.add_argument('--len_windows', '-l', required=True,
                        help='''
                        Size of the windows to be generated
                        ''', type=int, nargs=1)
    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument('--random', '-r', action='store_true',
                       help='''
                        Specifies program to create random, non-contiguous
                        alignment slices
                        ''')
    group.add_argument('--step_size', '-s', type=int, nargs='?',
                       const=25,
                       help='''
                       Defines the step size between contiguous windows -- 
                       default value is 25
                       ''')
    parser.add_argument('--num_windows', '-n', default=-1, nargs='?', 
                        const=500,
                        help='''
                        Defines number of windows to be created. --step_size
                        option overrides this value for contiguous windows. 
                        Default value for random windows is 500.
                        ''')
    parser.add_argument('--output_type', '-o', default='text_list',
                        const='fasta', type=str, nargs='?',
                        help='''
                        Defines how window slices are to be stored. Default
                        storage behavior is to store coordinates nucleotide
                        positions as a list of arrays in text file saved in 
                        {IN_FILE}-windows_list.txt. Option 'fasta' stores each
                        window as a separate fasta file in a new directory as 
                        {IN_FILE}-windowN_Stepk_Sizel.fas for contiguous
                        windows and {IN_FILE}-windowN_Sizel.fas for random
                        windows. NOTE: STORING AS SEPARATE FASTAS FOR RANDOM
                        WINDOWS IS NOT RECOMMENDED.
                        ''')
    args = vars(parser.parse_args(argv))
    in_file = args['in_file'][0]
    alignment = list(SeqIO.parse(in_file, "fasta"))
    aln_length = len(alignment[0].seq)
    len_windows = args['len_windows'][0]
    random = args['random']
    output = args['output_type']
    if random:
        num_windows = int(args['num_windows']) if int(args['num_windows']) > 0 else 500
        step = -1
    else:
        step = args['step_size']
        num_windows = math.floor((aln_length - len_windows)/step) + 1

    if random:
        windows = generate_random_windows(aln_length, len_windows,
                                          num_windows)
    else:
        windows = generate_contig_windows(len_windows, num_windows,
                                          step)
    if output == 'fasta':
        _store_as_fasta(alignment, windows, random, len_windows,
                        num_windows, step, in_file)
    else:
        _store_as_array_text(windows, random, len_windows,
                             num_windows, step, in_file)

    print('runtime complete')



if __name__ == "__main__":
    main(sys.argv[1:])
        
