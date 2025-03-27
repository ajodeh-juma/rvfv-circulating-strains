#!/usr/bin/env python3
"""
append addition columns to distinguish vaccine strains from other strains

usage:
python addStrainsTypes.py --metadata <file>  --prefix <str> --outDir <dir>
"""

import os
import sys
import logging
import argparse
import textwrap
from textwrap import dedent

import pandas as pd
import numpy as np
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
    required_group.add_argument('--alignment', required=True, type=str,
                                dest="alignment", metavar="<str>",
                                help="path to the alignment file FASTA format"
                                )
    required_group.add_argument('--lineage', required=True, type=str,
                                dest="lineage", metavar="<str>",
                                help="path to the file containing lineage assignment"
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
    required_group.add_argument('--prefix', required=True, metavar='<str>',
                                dest="prefix",
                                help="prefix of the output file(s)"
                                )
    parser.add_argument('--outDir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default=".")
    return parser


vaccine_strains = ['Smithburn', 'MP-12', '13', 'ZH-501-777', 'vacEGY', 'RVax-1', 'RSA/OBP/RVFVClone13/LAV',
                   'RSA/OBP/RVFV/Inactivated', 'RSA/OBP/RVFVSmithburn/LAV', 'rMP12-*NSs16/198', 'Clone 13',
                   '74HB59']

vaccine_accessions = ['DQ380193', 'OP146108', 'DQ380208', 'DQ380202', 'DQ380213', 'OP146105', 'DQ380212', 'OP146111']


def find_strain_type(x):
    if x in vaccine_strains:
    # if x in vaccine_accessions:
        return 'vaccine'
    else:
        return 'non-vaccine'


def strains_type(metadata, prefix, outdir):
    """
    add strains type

    :@param metadata:
    :@param prefix:
    :@param outdir:
    :return:
    """

    # create output directory if not exists
    outdir = os.path.abspath(outdir)
    mkdir(outdir)

    out_file = os.path.join(outdir, prefix) + '.csv'

    df = pd.read_csv(metadata, sep=",")
    cols = list(df.columns)
    logging.info("columns found {}".format(cols))
    df['strain_type'] = df['strain'].apply(find_strain_type)
    # df['strain_type'] = df['accession'].apply(find_strain_type)
    df.to_csv(out_file, sep=",", index=False)

    logging.info("output files: {}".format(out_file))
    return out_file


def merge_alignment_lineage_metadata(alignment, lineage, metadata, prefix, outdir):
    """

    :param alignment:
    :param lineage:
    :param metadata:
    :param prefix:
    :param outdir:
    :return:
    """

    # create output directory if not exists
    outdir = os.path.abspath(outdir)
    mkdir(outdir)

    outfile1 = os.path.join(outdir, prefix) + '.traits.txt'
    outfile2 = os.path.join(outdir, prefix) + '.fasta'

    # read metadata
    df_meta = pd.read_csv(metadata, sep=",")
    df_meta['year'] = df_meta.date.str.split('-').str[0]

    # read lineages
    df_lineage = pd.read_csv(lineage)

    # merge and sort by specified column
    meta = pd.merge(df_meta, df_lineage, left_on="accession", right_on="Query", how='inner').rename(columns={'Lineage': 'lineage'})
    meta['strain_type'] = meta['strain'].apply(find_strain_type)
    # meta['strain_type'] = meta['accession'].apply(find_strain_type)
    meta = meta[['accession', 'lineage', 'host', 'country', 'date', 'year', 'strain', 'strain_type']]
    meta['taxa'] = meta.accession.str.cat(meta.lineage, sep="|").str.cat(meta.host, sep="|").str.cat(meta.country,
                                                                                                     sep="|").str.cat(
        meta.strain_type, sep="|").str.cat(meta.year, sep="|")
    meta.insert(len(meta.columns) - 1, "taxa", meta.pop("taxa"))

    meta = meta.sort_values(by=['lineage', 'accession'], ascending=[True, True])
    # convert to dict
    d = meta.set_index('accession').T.to_dict('list')

    meta = meta[['taxa', 'accession', 'lineage', 'host', 'country', 'date', 'year', 'strain', 'strain_type']]

    # generator for multifasta file
    f_dic = dict()
    fasta_dict = dict(fasta_iterator(alignment))
    for k, v in fasta_dict.items():
        f_dic[k.split('|')[0]] = v

    # retain only accessions in the alignment/fasta file
    meta = meta[meta['accession'].isin(f_dic.keys())]

    with open(outfile2, 'w') as f_obj:
        for k, v in f_dic.items():
            if d.get(k) is None:
                continue
            else:
                f_obj.write(">" + str(d.get(k)[-1]) + "\n")
                f_obj.write(textwrap.fill(v, 80))
                f_obj.write("\n")

    print(meta.groupby(['lineage', 'strain_type']).count())
    meta.to_csv(outfile1, sep="\t", index=False)
    return outfile1, outfile2


def main():
    parser = parse_args()
    args = parser.parse_args()

    # check the arguments
    if args.metadata is not None and isinstance(args.metadata, str):
        logging.info("input metadata = {:^10}".format(os.path.basename(args.metadata), str))

    logging.info("appending additional strain information")
    # strains_type(metadata=args.metadata, prefix=args.prefix, outdir=args.outdir)
    merge_alignment_lineage_metadata(alignment=args.alignment,
                                     lineage=args.lineage,
                                     metadata=args.metadata,
                                     prefix=args.prefix,
                                     outdir=args.outdir)


if __name__ == '__main__':
    main()
