#!/usr/bin/env python3
"""
write each fasta sequence to a single file with the header as the filename

usage:
python write_each_sequence.py --fasta <fasta> --outDir .
"""

import os
import sys
import logging
import argparse
import textwrap
from textwrap import dedent

from utils import mkdir
from utils import fasta_iterator

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

    required_group.add_argument('--multi-fasta', required=True, type=str,
                                dest="multi_fasta", metavar="<str>",
                                help="path to the multi-fasta file with several sequences"
                                )
    parser.add_argument('--outDir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default=".")
    return parser


def write_fasta(multi_fasta, outdir):
    """

    @param multi_fasta:
    @param outdir:
    @return:
    """

    # create output directory if not exists
    outdir = os.path.abspath(outdir)
    mkdir(outdir)

    fasta_dict = dict(fasta_iterator(multi_fasta))
    for k, v in fasta_dict.items():
        outfile = os.path.join(outdir, k.split("|")[0] + '.fasta')
        if os.path.exists(outfile):
            continue
        else:
            logging.info("writing sequence accession {}".format(k))
            with open(outfile, 'w') as outfile:
                outfile.write(">" + k.split("|")[0] + "\n")
                outfile.write(textwrap.fill(v, 80))
                outfile.write("\n")
    return outdir


def main():
    parser = parse_args()
    args = parser.parse_args()

    write_fasta(multi_fasta=args.multi_fasta, outdir=args.outdir)


if __name__ == '__main__':
    main()
