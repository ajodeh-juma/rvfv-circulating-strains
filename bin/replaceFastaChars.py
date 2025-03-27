#!/usr/bin/env python3

import os
import sys
import logging
import argparse
import textwrap
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
from utils import mkdir

log_level = logging.DEBUG
logging.basicConfig(level=log_level, format='[%(asctime)s] - %(''levelname)s - [%(funcName)s] - %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S %p')


def parse_args():
    # arguments parser
    parser = argparse.ArgumentParser(
        prog="replaceFastaChars.py",
        argument_default=argparse.SUPPRESS,
        description=textwrap.dedent('''\
            replace characters in a multifasta file
            -------------------------------------------------------------------------

            '''),
    )
    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument('--fasta', required=True, type=str,
                                dest="fasta", metavar="<str>",
                                help="path to the multifasta FASTA format file having sequences to subset"
                                )
    required_group.add_argument('--outfile', required=False, metavar='<str>',
                                dest="outfile",
                                help="filename to the output file(s)"
                                )
    return parser


def replace_chars(fasta, outfile):
    """

    @param fasta:
    @param start:
    @param end:
    @param outfile:
    @return:
    """

    outdir = os.path.dirname(outfile)
    outdir = os.path.abspath(outdir)
    mkdir(outdir)

    with open(outfile, "wt") as out:
        for seq_record in SeqIO.parse(fasta, "fasta"):
            out.write(">" + str(seq_record.id) + "\n")
            out.write(str(seq_record.seq).replace('-', 'N') + "\n")

    logging.info("output files: fasta -> {}".format(outfile))
    return outfile


def main():
    parser = parse_args()
    args = parser.parse_args()

    # check the arguments
    if args.fasta is not None and isinstance(args.fasta, str):
        args.fasta = args.fasta
        logging.info("input multifasta = {:^10}".format(args.fasta, str))


    if args.outfile is not None and isinstance(args.outfile, str):
        args.outfile = args.outfile
        logging.info("output file = {:^10}".format(args.outfile, str))

    replace_chars(fasta=args.fasta, outfile=args.outfile)


if __name__ == '__main__':
    main()
