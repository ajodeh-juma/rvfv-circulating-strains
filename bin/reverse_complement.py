#!/usr/bin/env python3
"""
Reverse complement nucleotide sequence
"""

from __future__ import division
import os
import sys
import logging
import argparse
import textwrap

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

# import time
# import argparse
from os import path
# from datetime import datetime
import itertools
from textwrap import dedent
from utils import mkdir

log_level = logging.DEBUG
logging.basicConfig(level=log_level, format='[%(asctime)s] - %(''levelname)s - [%(funcName)s] - %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S %p')


def parse_args():
    # arguments parser
    parser = argparse.ArgumentParser(
        argument_default=argparse.SUPPRESS, prefix_chars='--',
        description=dedent(__doc__)
    )

    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument('--fasta', required=True, type=str,
                                dest="fasta", metavar="<str>",
                                help="path to the FASTA format file"
                                )
    required_group.add_argument('--prefix', required=True, metavar='<str>',
                                dest="prefix",
                                help="prefix of the output file(s)"
                                )
    parser.add_argument('--outDir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default=".")
    return parser


def reverse_complement(fasta, prefix, outdir):
    complement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                  'n': 'n', 'm': 'm', 'r': 'r', 'y': 'y',
                  'w': 'w', '-': '-', 's': 's', 'k': 'k',
                  'h': 'h'}

    fasta_dict = dict()
    revcomp_dict = dict()
    headers = 0
    for line in open(fasta):
        if line.startswith('>'):
            header = line.strip()
            headers += 1
            fasta_dict[header] = ''
        else:
            sequence = ''.join(line.strip())
            sequence = sequence.lower()
            fasta_dict[header] += sequence

      # create output directory if not exists
    outdir = os.path.abspath(outdir)
    mkdir(outdir)

    outfile = os.path.join(outdir, prefix) + '.reverse.complement.fasta'

    logging.info("Found {} sequences".format(len(fasta_dict.keys())))
    for seqid, seq in fasta_dict.items():
        revcomp = ''.join([complement[base] for base in seq[::-1]])
        fasta_dict[seqid] = revcomp

    with open(outfile, 'w') as f_obj:
        for k, v in fasta_dict.items():
            f_obj.write(k)
            f_obj.write('\n')
            if len(v) > 10:
                f_obj.write('\n'.join(v[i:i + 80] for i in range(0, len(v), 80)))
                f_obj.write('\n')
            else:
                f_obj.write(v)
                f_obj.write('\n')
    return outfile


def main():
    parser = parse_args()
    args = parser.parse_args()

    # check the arguments
    if args.fasta is not None and isinstance(args.fasta, str):
        args.fasta = args.fasta
        logging.info("input fasta = {:^10}".format(args.fasta, str))

    reverse_complement(fasta=args.fasta, prefix=args.prefix, outdir=args.outdir)


if __name__ == '__main__':
    main()
