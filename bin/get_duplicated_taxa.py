#!/usr/bin/env python3
"""
given output from seqkit rmdup duplicates ids file, the function extracts duplicated IDs
"""

import os
import sys
import logging
import argparse
import textwrap
import pycountry
from textwrap import dedent

import pandas as pd

from utils import mkdir
from utils import str2bool
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
    required_group.add_argument('--dup', required=True, type=str,
                                dest="dup", metavar="<str>",
                                help="path to the text file containing duplicated ids from seqkit rmdup command using "
                                     "-D option "
                                )
    required_group.add_argument('--outliers', required=True, type=str,
                                dest="outliers", metavar="<str>",
                                help="path to the text file containing outlier sequence identifiers, one per line"
                                )
    required_group.add_argument('--prefix', required=True, metavar='<str>',
                                dest="prefix",
                                help="prefix of the output file(s)"
                                )
    parser.add_argument('--remove-duplicates', metavar="<boolean>", type=str2bool, default=False,
                        dest='remove_duplicates',
                        help=dedent('''remove duplicated sequences from the dataset''')
                        )
    parser.add_argument('--outDir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default=".")
    return parser


def list_dups(dup, prefix, outliers, remove_duplicates, outdir):
    """

    :param dup:
    :param prefix:
    :param outdir:
    :return:
    """

    # create output directory if not exists
    outdir = os.path.abspath(outdir)
    mkdir(outdir)

    out_file = os.path.join(outdir, prefix) + '_duplicated_taxa.txt'

    # read dataframe
    df = pd.read_csv(dup, sep="\t", header=None, names=['number', 'taxa'])
    df['duplicates'] = df.taxa.str.split(',').str[1:]
    taxa = [x for xs in list(df.duplicates) for x in xs]

    if remove_duplicates:
        taxa = taxa
    else:
        taxa = []

    logging.info("excluding duplicates: {}".format(len(taxa)))


    # read outliers
    outliers = os.path.abspath(outliers)
    with open(outliers) as fh:
        newlist = [line.rstrip() for line in fh.readlines()]
    
    # concatenate both lists
    all_taxa = taxa + newlist

    logging.info("duplicated and outlier taxa {}".format(len(all_taxa)))
    with open(out_file, 'w') as f_obj:
        for t in all_taxa:
            f_obj.write(t.strip()+'\n')
    return out_file


def main():
    parser = parse_args()
    args = parser.parse_args()

    # check the arguments
    if args.dup is not None and isinstance(args.dup, str):
        args.dup = args.dup
        logging.info("input dup = {:^10}".format(args.dup, str))

    logging.info("listing duplicated taxa")
    list_dups(dup=args.dup, prefix=args.prefix, outliers=args.outliers, remove_duplicates=args.remove_duplicates, outdir=args.outdir)


if __name__ == '__main__':
    main()
