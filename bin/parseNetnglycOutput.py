#!/usr/bin/env python3
"""
filter snps based on threshold of occurence in sequences
"""

import os
import csv
import sys
import logging
import argparse
import textwrap
from textwrap import dedent
import pandas as pd

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
    required_group.add_argument('--input',
                                required=True,
                                type=str,
                                dest="input",
                                metavar="<str>",
                                help="output file from netNglyc analysis"
                                )
    required_group.add_argument('--prefix',
                                required=True,
                                metavar='<str>',
                                dest="prefix",
                                help="prefix of the output file(s)"
                                )
    parser.add_argument('--outdir',
                        metavar='<DIR>',
                        dest="outdir",
                        default=".",
                        help="path to the output directory")
    return parser


def parse_netglyc(input, prefix, outdir):
    """

    :param netglyc:
    :param prefix:
    :param outdir:
    :return:
    """

    c = 0
    b = False
    record = []
    records = []
    with open(input) as f_input:
        for line in f_input:
            if 'agreement result' in line:
                b = True
                continue
            elif line.startswith('Graphics in PostScript'):
                b = False
                continue
            elif b:
                if line.startswith('-'):
                    continue
                else:
                    records += [line.strip().split()]
    # write to file
    outdir = os.path.abspath(outdir)
    mkdir(outdir)

    outfile = os.path.join(outdir, prefix + '.csv')

    # header = ['SeqName', 'Position', 'Potential', 'Jury', 'Agreement',
    # 'Result']
    header = ['Accession', 'Host', 'Country', 'Date', 'Position', 'Potential', 'Jury', 'Agreement', 'Result']
    d = dict()
    with open(outfile, 'w', encoding='UTF8', newline='') as f_obj:
        writer = csv.writer(f_obj)
        writer.writerow(header)
        for i in range(len(records)):
            if len(records[i]) == 0 or records[i][5] == '-':
                continue
            else:
                seqname = records[i][0].split("_")
                if 'NC' in seqname:
                    accession = seqname[0]+'_'+seqname[1]
                    host = seqname[2]
                    country = seqname[3]
                    date = seqname[4]
                else:
                    accession = seqname[0]
                    host = seqname[1]
                    country = seqname[2]
                    date = seqname[3]
                new_records = [accession, host, country, date]+records[i][1:6]
                # writer.writerow(records[i][0:6])
                writer.writerow(new_records)
    return outfile


def main():
    parser = parse_args()
    args = parser.parse_args()

    # check the arguments
    if args.input is not None and isinstance(args.input, str):
        logging.info("input file = {:^10}".format(os.path.basename(args.input), str))

    logging.info("parsing netNglyc {}".format(args.input))

    parse_netglyc(input=args.input, prefix=args.prefix, outdir=args.outdir)


if __name__ == '__main__':
    main()