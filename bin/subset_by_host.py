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
    # required_group.add_argument('--host', required=True, type=str,
    #                             dest="host", metavar="<str>",
    #                             choices=['bat', 'buffalo', 'cow', 'goat', 'human', 'mosquito', 'sheep', 'springbok'],
    #                             help="host name (must be in the metadata csv and fasta file)"
    #                             )
    required_group.add_argument('--prefix', required=True, metavar='<str>',
                                dest="prefix",
                                help="prefix of the output file(s)"
                                )
    parser.add_argument('--outDir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default=".")
    return parser


def subset_data(fasta, csv, prefix, outdir):
    """
    subset sequences by host

    :param fasta:
    :param csv:
    # :param host:
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
    df = pd.read_csv(csv, sep="\t")

    # iterate through the input fasta file
    fasta_dict = dict(fasta_iterator(fasta))

    # subset host
    hosts = list(set(list(df.host)))
    for h in hosts:
        fname = '_'.join(h.split())
        out1 = os.path.join(outdir, prefix + '.' + fname) + '.csv'
        out2 = os.path.join(outdir, prefix + '.' + fname) + '.fasta'
        subdf = df[df['host'].isin([h])]
        subdf.to_csv(out1, sep=",", index=False)
        sub_d = subdf.set_index('taxa').T.to_dict('list')
        logging.info("{} host found {} sequences".format(h, len(subdf)))
        with open(out2, 'w') as outfile:
            for k, v in sub_d.items():
                outfile.write(">" + k + "\n")
                outfile.write(textwrap.fill(fasta_dict.get(k), 80))
                outfile.write("\n")

    # subset_df = df[df['host'].isin([host])]
    # subset_d = subset_df.set_index('taxa').T.to_dict('list')
    # subset_df.to_csv(out_csv, sep=",", index=False)

    # logging.info("found {} sequences".format(len(subset_df)))

    # # get and write to file
    # with open(out_fasta, 'w') as outfile:
    #     for k, v in subset_d.items():
    #         outfile.write(">" + k + "\n")
    #         outfile.write(textwrap.fill(fasta_dict.get(k), 80))
    #         outfile.write("\n")
    # return out_csv, out_fasta


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

    subset_data(fasta=args.fasta, csv=args.csv, prefix=args.prefix, outdir=args.outdir)


if __name__ == '__main__':
    main()
