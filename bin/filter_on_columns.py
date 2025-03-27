#!/usr/bin/env python3
"""
filter out sequences based on metadata columns
"""

import os
import logging
import argparse
import textwrap
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
    required_group.add_argument('--columns', required=True, nargs='+',
                                dest="columns", metavar="<str>",
                                help="column names (separated by space) to exclude where thera are NAs values (i.e "
                                     "date, host, country, location. The date column should be in ISO8601 format")
    required_group.add_argument('--prefix', required=True, metavar='<str>',
                                dest="prefix",
                                help="prefix of the output file(s)"
                                )
    parser.add_argument('--outDir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default=".")
    return parser


def subset_data(fasta, csv, columns, prefix, outdir):
    """
    filter out sequences without locations

    :param fasta:
    :param csv:
    :param columns:
    :param prefix:
    :param outdir:
    :return:
    """

    # create output directory if not exists
    outdir = os.path.abspath(outdir)
    mkdir(outdir)

    out_csv = os.path.join(outdir, prefix) + '.csv'
    out_fasta = os.path.join(outdir, prefix) + '.fasta'

    df = pd.read_csv(csv, sep=",")
    logging.info("found {} sequences".format(len(df)))
    # filter out sequences on given columns
    columns = columns.split()
    df = df.dropna(subset=columns, how='any')
    # identify duplicates and drop
    duplicate = df[df.duplicated('accession')]
    df = df.sort_values(by=['accession', 'length'])
    df = df.drop_duplicates(subset=['accession'], keep='last')
    
    # new dataframe with column taxa similar to fasta headers
    # accession,organism,length,strain,country,location,host,date
    new_df = df.copy()
    new_df['taxa'] = df[['accession', 'host', 'country', 'location', 'date']].fillna('').agg('|'.join, axis=1)
    # new_df['taxa'] = new_df.accession.str.cat([new_df.host, new_df.country, new_df.location, new_df.date], sep="|")
    new_df.reset_index(drop=True, inplace=True)
    d = new_df.set_index('taxa').T.to_dict('list')
    logging.info("subset {} sequences".format(len(d)))

    # iterate through the input fasta file
    fasta_dict = dict(fasta_iterator(fasta))
    f = []
    with open(out_fasta, 'w') as outfile:
        for k, v in d.items():
            if fasta_dict.get(k) is None:
                # print(k, fasta_dict.get(k.rsplit("|")[0]))
                continue
            else:
                f.append(k.rsplit("|")[0])
                outfile.write(">" + k + "\n")
                outfile.write(textwrap.fill(fasta_dict.get(k), 80))
                outfile.write("\n")

    df = new_df[new_df['accession'].isin(f)]
    df = df.drop(['taxa'], axis=1)
    df.to_csv(out_csv, sep=",", index=False)
    return out_csv, out_fasta


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

    if args.columns is not None:
        columns = ' '.join(str(k) for k in args.columns)
        logging.info("selected column = {:^10}".format(columns))

    logging.info("Sub-setting sequences on column: {}".format(args.columns))
    subset_data(fasta=args.fasta, csv=args.csv, columns=columns, prefix=args.prefix, outdir=args.outdir)


if __name__ == '__main__':
    main()
