#!/usr/bin/env python3
"""
filter out outlier sequences detected by tempest
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
    required_group.add_argument('--meta', required=True, type=str,
                                dest="meta", metavar="<str>",
                                help="path to the file containing the FASTA sequences metadata"
                                )
    required_group.add_argument('--sep', required=True, type=str,
                                dest="sep", metavar="<str>", choices=['comma', 'tab'],
                                help="separator used in the metadata file"
                                )
    required_group.add_argument('--taxa', required=True, type=str,
                                dest="taxa", metavar="<str>",
                                help="path to the text file containing taxa names per line"
                                )
    required_group.add_argument('--prefix', required=True, metavar='<str>',
                                dest="prefix",
                                help="prefix of the output file(s)"
                                )
    parser.add_argument('--outDir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default=".")
    return parser


def filter_taxa(fasta, meta, sep, taxa, prefix, outdir):
    """
    filter out outlier sequences detected by tempest

    :@param fasta:
    :@param meta:
    :@param sep:
    :@param taxa:
    :@param prefix:
    :@param outdir:
    :return:
    """

    # create output directory if not exists
    global ext
    outdir = os.path.abspath(outdir)
    mkdir(outdir)

    if sep == 'comma':
        sep = ','
        ext = '.csv'
    elif sep == 'tab':
        sep = '\t'
        ext = '.txt'
    else:
        logging.info("unknown separator")

    out_file = os.path.join(outdir, prefix) + '_filtered' + ext
    out_fasta = os.path.join(outdir, prefix) + '_filtered.fasta'

    lines = []
    with open(taxa) as f:
        for line in f:
            lines.append(line.strip())
    lines = list(set(lines))

    # read dataframe
    df = pd.read_csv(meta, sep=sep)

    # filter out sequences
    df1 = df[~df['taxa'].isin(lines)]
    cols = list(df.columns)
    if 'country' in cols:
        df2 = df1.sort_values(by=['country'], ascending=True)
    elif 'host' in cols:
        df2 = df1.sort_values(by=['host'], ascending=True)
    else:
        df2 = df1.sort_values(by=['taxa'], ascending=True)

    # iterate through the input fasta file
    fasta_dict = dict(fasta_iterator(fasta))

    # write output
    with open(out_fasta, 'w') as outfile:
        for acc in list(df2.taxa):
            # if fasta_dict.get(acc) is None:
            #     continue
            # else:
            outfile.write(">" + str(acc) + "\n")
            outfile.write(fasta_dict.get(acc))
            # outfile.write(textwrap.fill(fasta_dict.get(acc), 80))
            outfile.write("\n")

    ls1 = fasta_dict.keys()
    dups = list(set(list(ls1)) ^ set(lines))
    outliers = list(set(list(df2.taxa)) ^ set(ls1))

    logging.info("total taxa found: {}".format(len(fasta_dict)))
    # logging.info("duplicated taxa: {}".format(len(lines)))
    logging.info("outlier taxa: {}".format(len(outliers)))
    logging.info("new taxa: {}".format(len(df2)))

    df2.to_csv(out_file, sep=sep, index=False)
    logging.info("output files:\nmetadata -> {}\nfasta -> {}".format(out_file, out_fasta))
    return out_file, out_fasta


def main():
    parser = parse_args()
    args = parser.parse_args()

    # check the arguments
    if args.fasta is not None and isinstance(args.fasta, str):
        args.fasta = args.fasta
        logging.info("input fasta = {:^10}".format(args.fasta, str))

    if args.meta is not None and isinstance(args.meta, str):
        args.meta = args.meta
        logging.info("input metadata file = {:^10}".format(args.meta, str))

    if args.taxa is not None and isinstance(args.taxa, str):
        args.taxa = args.taxa
        logging.info("taxa to exclude file = {:^10}".format(args.taxa, str))
    # elif not isinstance(args.country, str):
    #     args.country = repr(args.country)
    # else:
    #     print("please specify the country in string format")
    #     sys.exit(2)

    logging.info("filtering outlier taxa")
    filter_taxa(fasta=args.fasta, meta=args.meta, sep=args.sep, taxa=args.taxa, prefix=args.prefix, outdir=args.outdir)


if __name__ == '__main__':
    main()
