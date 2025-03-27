#!/usr/bin/env python3
"""
rename headers and generate a corresponding tab-delimited file with new headers

usage:
python reformat_headers.py --fasta <fasta> --csv <csv> --headers <accession year> --prefix <acc.year> --outDir .
"""

import os
import sys
import logging
import argparse
import textwrap
from textwrap import dedent

import pandas as pd
import pycountry

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
    required_group.add_argument('--fasta', required=True, type=str,
                                dest="fasta", metavar="<str>",
                                help="path to the FASTA format file"
                                )
    required_group.add_argument('--metadata', required=True, type=str,
                                dest="metadata", metavar="<str>",
                                help="path to the corresponding metadata file containing the columns [ "
                                     "organism,length,strain,country,location,host,date] for comma-separated format "
                                     "file or [lat long] for coordinates tab delimited file. All these file formats "
                                     "must contain in the FIRST column containing sequence identifiers similar to the "
                                     "ones in the "
                                     "fasta file"
                                )
    required_group.add_argument('--sep', required=True, type=str,
                                dest="sep", metavar="<str>", choices=['comma', 'tab'],
                                help="separator used in the metadata file"
                                )
    required_group.add_argument('--headers', required=True, nargs='+',
                                dest="headers", metavar="<str>",
                                choices=['year', 'date', 'host', 'country', 'location'],
                                help="space-separated list of the type of header details to include (i.e "
                                     "date, year, host, country. The date column should be in ISO8601 format")
    required_group.add_argument('--prefix', required=True, metavar='<str>',
                                dest="prefix",
                                help="prefix of the output file(s)"
                                )
    parser.add_argument('--coordinates', type=str2bool, metavar="<bool>",
                        default=True, const=True, nargs='?',
                        dest='coordinates',
                        help="specify if the metadata file contains coordinates with columns 'lat' and 'long'"
                        )
    parser.add_argument('--outDir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default=".")
    return parser


def rename_headers(fasta, metadata, sep, headers, prefix, coordinates, outdir):
    """
    rename headers and generate a corresponding tab-delimited file with new headers

    :param fasta:
    :param metdata:
    :@param headers:
    :param prefix:
    :@param coordinates:
    :param outdir:
    :return:
    """

    # check separator
    if sep == 'comma':
        sep = ','
        ext = '.csv'
    elif sep == 'tab':
        sep = '\t'
        ext = '.txt'
    else:
        logging.info("unknown separator")



    # countries
    # alpha2ctry = dict()
    # t = list(pycountry.countries)
    # for country in t:
    #     alpha2ctry[country.alpha_3] = country.name.split(',')[0]

    # create output directory if not exists
    outdir = os.path.abspath(outdir)
    mkdir(outdir)

    out_file = os.path.join(outdir, prefix) + '.txt'
    out_fasta = os.path.join(outdir, prefix) + '.fasta'

    df = pd.read_csv(metadata, sep=sep)
    cols = list(df.columns)
    logging.info("columns found {}".format(cols))
    if coordinates:
        expected = ['lat', 'long']
        check = all(item in cols for item in expected)
        if check:
            # subset dataframe to exclude first column (accession)
            old_df = df[list(df.columns)[1:]]
            # create columns
            df['accession'] = df.traits.str.split("|").str[0]
            df['host'] = df.traits.str.split("|").str[1]
            df['country'] = df.traits.str.split("|").str[2]
            df['location'] = df.traits.str.split("|").str[3]
            df['date'] = df.traits.str.split("|").str[4]
            df['year'] = df.traits.str.split("|").str[4].str.split('-').str[0]
            df['taxa'] = df.accession.str.cat([df[k] for k in headers], sep="|")
            # df['date'] = pd.to_datetime(df['date'], format="%d-%m-%Y") 
            d = df.set_index('accession').T.to_dict('list')
            headers.insert(0, 'taxa')
            df = df[['taxa', 'lat', 'long']]
            df.to_csv(out_file, sep="\t", index=False)
    else:
        expected = ['date', 'host', 'country', 'location']
        check = all(item in cols for item in expected)
        if check:
            # sort and drop duplicates
            df = df.sort_values(by=['accession', 'length'])           
            df = df.drop_duplicates(subset=['accession'], keep='last')
            df['year'] = df.date.str.split('-').str[0]
            df['taxa'] = df.accession.str.cat([df[k] for k in headers], sep="|")
            df.insert(len(df.columns)-1, "taxa", df.pop("taxa"))

            # convert to dict
            d = df.set_index('accession').T.to_dict('list')
            headers.insert(0, 'taxa')
            df = df[headers[:-1]]
            df.to_csv(out_file, sep="\t", index=False)

    # generator for multifasta file
    f_dic = dict()
    fasta_dict = dict(fasta_iterator(fasta))
    for k, v in fasta_dict.items():
        f_dic[k.split('|')[0]] = v

    with open(out_fasta, 'w') as outfile:
        for k, v in f_dic.items():
            if d.get(k) is None:
                continue
            else:
                outfile.write(">" + str(d.get(k)[-1]) + "\n")
                outfile.write(textwrap.fill(v, 80))
                outfile.write("\n")

    logging.info("output files:\ncsv -> {}\nfasta -> {}".format(out_file, out_fasta))
    return out_file, out_fasta


def main():
    parser = parse_args()
    args = parser.parse_args()

    # check the arguments
    if args.fasta is not None and isinstance(args.fasta, str):
        logging.info("input fasta = {:^10}".format(args.fasta, str))
    if args.metadata is not None and isinstance(args.metadata, str):
        logging.info("input metadata = {:^10}".format(args.metadata, str))
    if args.headers is not None:
        headers = ' '.join(str(k) for k in args.headers)
        logging.info("selected headers = {:^10}".format(headers))

    logging.info("reformating headers")
    rename_headers(
        fasta=args.fasta, metadata=args.metadata, headers=args.headers, sep=args.sep,
        prefix=args.prefix, coordinates=args.coordinates, outdir=args.outdir
    )


if __name__ == '__main__':
    main()
