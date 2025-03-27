#!/usr/bin/env python3
"""
subset sequences using a start and end position
"""

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
        prog="subsetByCoordinates.py",
        argument_default=argparse.SUPPRESS,
        description=textwrap.dedent('''\
            extract sequences by specified coordinates (start and end positions from a multifasta)
            -------------------------------------------------------------------------
            start integer
            end integer
            '''),
    )
    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument('--fasta', required=True, type=str,
                                dest="fasta", metavar="<str>",
                                help="path to the multifasta FASTA format file having sequences to subset"
                                )
    required_group.add_argument('--start', required=True, type=int,
                                dest="start", metavar="<int>",
                                help="integer specifying the start position (1-based indexing)"
                                )
    required_group.add_argument('--end', required=True, type=int,
                                dest="end", metavar="<int>",
                                help="integer specifying the end position (1-based indexing)"
                                )
    required_group.add_argument('--outfile', required=False, metavar='<str>',
                                dest="outfile",
                                help="filename to the output file(s)"
                                )
    return parser


def subset_by_position(fasta, start, end, outfile):
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
            out.write(str(seq_record.seq[int(start):int(end)]) + "\n")

    logging.info("output files: fasta -> {}".format(outfile))
    return outfile


def main():
    parser = parse_args()
    args = parser.parse_args()

    # check the arguments
    if args.fasta is not None and isinstance(args.fasta, str):
        args.fasta = args.fasta
        logging.info("input multifasta = {:^10}".format(args.fasta, str))

    if args.start is not None and isinstance(args.start, int):
        args.start = args.start
        logging.info("start position = {:^10}".format(args.start, int))

    if args.end is not None and isinstance(args.end, int):
        args.end = args.end
        logging.info("end position = {:^10}".format(args.end, int))

    if args.start >= args.end:
        logging.info("start position {} cannot be greater or equal to end position {}".format(args.start, args.end))
    if args.end <= args.start:
        logging.info("end position {} cannot be less or equal to start position {}".format(args.end, args.start))

    if args.outfile is not None and isinstance(args.outfile, str):
        args.outfile = args.outfile
        logging.info("output file = {:^10}".format(args.outfile, str))

    logging.info("Sub-setting sequences on coordinates/positions {} and {}".format(args.start, args.end))
    subset_by_position(fasta=args.fasta, start=args.start, end=args.end, outfile=args.outfile)


if __name__ == '__main__':
    main()
