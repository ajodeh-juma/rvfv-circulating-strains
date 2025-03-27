#!/usr/bin/env python3
"""
subset sequences by host
"""

import os
import sys
import logging
import argparse
import textwrap
# import pycountry
from textwrap import dedent

import pandas as pd

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
    required_group.add_argument('--fasta', required=True, type=str,
                                dest="fasta", metavar="<str>",
                                help="path to the FASTA format file"
                                )
    required_group.add_argument('--csv', required=True, type=str,
                                dest="csv", metavar="<str>",
                                help="path to the CSV file containing the FASTA sequences metadata"
                                )
    required_group.add_argument('--years', required=True, type=int,
                                dest="years", metavar="<str>", nargs="+",
                                help="space separated list of years to subset data and extract"
                                )
    required_group.add_argument('--prefix', required=True, metavar='<str>',
                                dest="prefix",
                                help="prefix of the output file(s)"
                                )
    parser.add_argument('--outDir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default=".")
    return parser


def subset_data(fasta, csv, years, prefix, outdir):
    """
    subset sequences by host

    :param fasta:
    :param csv:
    :param years:
    :param prefix:
    :param outdir:
    :return:
    """

    # create output directory if not exists
    outdir = os.path.abspath(outdir)
    mkdir(outdir)

    out_csv = os.path.join(outdir, prefix) + '.csv'
    out_fasta = os.path.join(outdir, prefix) + '.fasta'

    # read dataframe
    df = pd.read_csv(csv, sep=",")
    df['date'] = pd.to_datetime(df['date'], format="%Y-%m-%d")
    df = df.fillna("")
    df['taxa'] = df["accession"].str.cat(df[["host", "country", "location", "date"]].astype(str), sep="|")
    df['year'] = df.date.dt.year

    # iterate through the input fasta file
    fasta_dict = dict(fasta_iterator(fasta))

    df = df[df['year'].isin(years)]

    sub_d = df.set_index('taxa').T.to_dict('list')
    logging.info("found {} sequences".format(len(df)))

    # write fasta
    with open(out_fasta, 'w') as outfile:
        for k, v in sub_d.items():
            outfile.write(">" + k + "\n")
            outfile.write(textwrap.fill(fasta_dict.get(k), 80))
            outfile.write("\n")
        
    # write csv metadata
    subdf = df.drop(['taxa', 'year'], axis=1)
    subdf.to_csv(out_csv, sep=",", index=False)


def main():
    parser = parse_args()
    args = parser.parse_args()

    # check the arguments
    if args.fasta is not None and isinstance(args.fasta, str):
        args.fasta = args.fasta
        logging.info("input fasta = {:^10}".format(args.fasta, str))

    if args.csv is not None and isinstance(args.csv, str):
        args.csv = args.csv
        logging.info("input csv = {:^10}".format(args.csv, str))

    # if args.host is not None:
    #     args.host = args.host
    #     logging.info("selected host = {:^10}".format(args.host))

    subset_data(fasta=args.fasta, csv=args.csv, years=args.years, prefix=args.prefix, outdir=args.outdir)


if __name__ == '__main__':
    main()
