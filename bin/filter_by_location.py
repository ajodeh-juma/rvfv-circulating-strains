#!/usr/bin/env python3
"""
filter out sequences without locations
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
    required_group.add_argument('--csv', required=True, type=str,
                                dest="csv", metavar="<str>",
                                help="path to the CSV file containing the FASTA sequences metadata"
                                )
    required_group.add_argument('--prefix', required=True, metavar='<str>',
                                dest="prefix",
                                help="prefix of the output file(s)"
                                )
    parser.add_argument('--outDir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default=".")
    return parser


def subset_data(fasta, csv, prefix, outdir):
    """
    filter out sequences without locations

    :param fasta:
    :param csv:
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

    # filter out sequences without locations
    # df = df[df['location'].str.len() > 0]
    # print(df.query('location != location'))
    print(df.loc[df['location'].isnull()])
    
    df = df[df['location'].notnull()]

    # identify duplicates and drop
    duplicate = df[df.duplicated('accession')]
    df = df.sort_values(by=['accession', 'length'])
    df = df.drop_duplicates(subset=['accession'], keep='last')


    # print(len(df))
    # accessions = list(df.accession)
    logging.info("found {} sequences".format(len(df)))
    # new dataframe with column taxa similar to fasta headers
    # accession,organism,length,strain,country,location,host,date
    new_df = df.copy()
    # >MSA|human|KEN|0-0-2007
    new_df['taxa'] = new_df.accession.str.cat([new_df.host, new_df.country, new_df.location, new_df.date], sep="|")
    new_df.reset_index(drop=True, inplace=True)
    d = new_df.set_index('taxa').T.to_dict('list')

    print(len(d))

    # iterate through the input fasta file
    fasta_dict = dict(fasta_iterator(fasta))
    f = []
    with open(out_fasta, 'w') as outfile:
        for k, v in d.items():
            if fasta_dict.get(k) is None:
                continue
            else:
                f.append(k)
                outfile.write(">" + k + "\n")
                outfile.write(textwrap.fill(fasta_dict.get(k), 80))
                outfile.write("\n")

    # print(len(d.keys()), len(sorted(f)))
    # for i, (j, k) in enumerate(zip(sorted(d.keys()), sorted(f))):
    #     if j != k:
    #         print(i, j, k)
    # #
    # # print(len(first))
    #
    # l = [k for k in d.keys()]
    df = new_df[new_df['taxa'].isin(f)]
    df.drop(['taxa'], axis=1)
    df.to_csv(out_csv, sep=",", index=False)

    # for k, v in fasta_dict.items():
    #     if k.split('|')[0] in accessions:
    #         f.append(k)
    #         first.append(k.split('|')[0])
    #
    # with open(out_fasta, 'w') as outfile:
    #     for acc in f:
    #         outfile.write(">" + acc + "\n")
    #         outfile.write(textwrap.fill(fasta_dict.get(acc), 80))
    #         outfile.write("\n")
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

    # if args.country is not None and isinstance(args.country, str):
    #     args.country = args.country
    #     print("country = {:^10}".format(args.country, str))
    # elif not isinstance(args.country, str):
    #     args.country = repr(args.country)
    # else:
    #     print("please specify the country in string format")
    #     sys.exit(2)

    logging.info("Sub-setting sequences on locations")
    subset_data(fasta=args.fasta, csv=args.csv, prefix=args.prefix, outdir=args.outdir)


if __name__ == '__main__':
    main()
