# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 21:44:53 2020

@author: Michael
"""

import argparse
import re
import sys
import yaml

from Bio import SeqIO


def _coord_type(arg_value):
    pat = re.compile(r"(\d+)\:(\d+)")
    match = pat.match(arg_value)
    if match:
        if (int(match.group(1)) < int(match.group(2))) and \
            int(match.group(1)) > 0 and \
            int(match.group(2)) > 0:
            return arg_value
        else:
            raise argparse.ArgumentTypeError(
                "Invalid input: {0}".format(arg_value)
            )

    else:
        raise argparse.ArgumentTypeError
    
def parse_cli(argv):
    parser = argparse.ArgumentParser(
        prog='SARS_subgenome_generator.py',
        description="""
        Grabs desired subgenomic sequences from
        supplied input file. Coordinates are gotten from default
        or supplied yaml file with named region or :-deliminted
        position list. Reference sequence is assumed to be sequence 1.
        """,
        epilog="""
        Any questions or bugs may be directed to michaeljbale93@gmail.com
        """,
        add_help=True
    )
    
    subparsers = parser.add_subparsers(help='sub-command help')
    yaml_parser = subparsers.add_parser('with_yaml', help="with_yaml help")
    yaml_parser.add_argument(
        '-y',
        '--yaml',
        required=True,
        help="""
        Input yaml file with named genomic regions
        """
    )
    yaml_parser.add_argument(
        '-r',
        '-n',
        '--region',
        default='all',
        help="""
        Region name to pull from supplied yaml file. Default is to parse
        through all regions in yaml.
        """
    )
    yaml_parser.add_argument(
        '-i',
        '--in_fasta',
        required=True,
        help="""
        Input fasta file with sequences to slice. 
        Format is assumed to be fasta.
        """
    )
    yaml_parser.add_argument(
        '-p',
        '--prefix',
        required=True,
        help="""
        Prefix name for output sequence file. 
        Output filenames will be generated as ${OUTPUT}_${REGION}.fas
        """
    )
    
    nt_coord_parser = subparsers.add_parser('nt_coord', help='nt_coord help')
    
    nt_coord_parser.add_argument(
        '-c',
        '--coords',
        required=True,
        type=_coord_type,
        help="""
        Coordinate positions to get from supplied sequence file. Must be in
        form XXX:YYY where XXX < YYY and both are greater than 0.
        """
    )    
    
    nt_coord_parser.add_argument(
        '-f',
        '--fasta',
        required=True,
        help="""
        Input fasta file with sequences to slice. 
        Format is assumed to be fasta.
        """
    )
    nt_coord_parser.add_argument(
        '-o',
        '--output',
        required=True,
        help="""
        Prefix name for output sequence file. 
        Output filenames will be generated as ${OUTPUT}_${REGION}.fas
        """
    )
    
    return parser.parse_args(argv)

def gather_subgenomes_by_yaml(args):
    with open(args['yaml'], 'r') as yml:
        def_regions = yaml.load(yml, Loader = yaml.FullLoader)
        
    if args['region'] == 'all':
        regions = [k for k in def_regions.keys()]
    elif args['region'] in def_regions.keys():
        regions = [args['region']]
    else:
        raise KeyError
    
    for region in regions:
        coords = "{0}:{1}".format(
            def_regions[region]['start'],
            def_regions[region]['end']
        )
        tmp_args = [
            'nt_coord',
            '-f', args['in_fasta'],
            '-c', coords,
            '-o', "{pref}_{name}".format(pref = args['prefix'], name=region)
        ]
        tmp_parsed = vars(parse_cli(tmp_args))
        gather_subgenome_by_coords(tmp_parsed, use_coords=False)
        
    
    
            
    
        

def gather_subgenome_by_coords(args, use_coords = True):
    seqs = list(SeqIO.parse(args['fasta'], "fasta"))
    glb_ref = seqs[0]
    nt_start = int(re.match((r"(\d+)\:(\d+)"), args['coords']).group(1)) - 1
    
    new_start = [
        m.start() for m in re.finditer(
            r"[^-]",
            str(glb_ref.seq)
        )
    ][nt_start]
    nt_end = int(re.match((r"(\d+)\:(\d+)"), args['coords']).group(2)) - 1
    
    new_end = [
        m.start() for m in re.finditer(
            r"[^-]",
            str(glb_ref.seq)
        )
    ][nt_end] + 1
    new_seqs = [seq[new_start:new_end] for seq in seqs]
    
    if use_coords:
        fname = "{}_coords_{}_{}.fas".format(args['output'], nt_start, nt_end)
    else:
        fname = args['output'] + ".fas"
    SeqIO.write(new_seqs, fname, "fasta")
    
    


def main(argv):
    args = vars(parse_cli(argv))

    if 'yaml' in args.keys():
        gather_subgenomes_by_yaml(args)
    else:
        gather_subgenome_by_coords(args)
        
if __name__ == "__main__":
    main(sys.argv[1:])