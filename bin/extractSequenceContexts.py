#!/usr/bin/env python3
"""
Extract flanking sequences for APOBEC sequence contexts
"""

import os
import re
import csv
import logging
import argparse
import textwrap
import pandas as pd
from textwrap import dedent
from itertools import groupby
from collections import defaultdict
from collections import Counter

from utils import str2bool
from utils import mkdir

log_level = logging.DEBUG
logging.basicConfig(level=log_level, format='[%(asctime)s] - %(''levelname)s - [%(funcName)s] - %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S %p')


def parse_args():
    # arguments parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        argument_default=argparse.SUPPRESS, prefix_chars='--',
        description=dedent('''get SNP records''')
    )
    # helpstr = """python3 fetch_segments.py [options]"""

    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument('--alignment', required=True, type=str,
                                dest="alignment", metavar="<str>",
                                help="path to the alignment file in fasta format"
                                )
    required_group.add_argument('--snps', metavar="<str>", dest="snps", default='',
                                help="a comma-separated file containing SNPs, must have the columns labelled 'snps', "
                                     "'positions', 'accession'")
    required_group.add_argument('--columns', required=True, nargs='+',
                                dest="columns", metavar="<str>",
                                choices=['country', 'date', 'host', 'lineage', 'strain_type', 'year'],
                                help="space-separated list of the type of header details to include (i.e "
                                     "date, year, host, country. The date column should be in ISO8601 format")
    required_group.add_argument('--mutation', required=True, metavar='<str>',
                                dest="mutation", choices=['CT', 'GA'],
                                help="mutation type to extract sequences"
                                )
    parser.add_argument('--outdir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default=".")
    return parser


def fasta_iterator(fasta_name):
    """
    Given a fasta file. yield tuples of header, sequence
    Author: Brent Pedersen
    https://www.biostars.org/p/710/
    """
    fa_fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fa_fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        header = next(header)[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in next(faiter))

        try:
            yield header, seq
        except StopIteration:
            return

        # yield header, seq
    fa_fh.close()


def extract_sequences(alignment, snps, mutation, columns, outdir):
    """

    :param alignment:
    :param snps:
    :param mutation:
    :param columns:
    :param outdir:
    :return:
    """

    df = pd.read_csv(snps)
    df['SNP'] = df['snps'].apply(lambda x: ''.join(re.findall(r"([A-Z]+)", str(x))))
    # df = df.assign(SNP=df["SNP"].apply(lambda l: "".join(l)))
    df1 = df[(df['SNP'] == mutation)]
    df1 = df[(df['strain_type'] == 'non-vaccine')]

    cols = ['accession', 'positions', 'SNP']

    ord_dict = df1[cols].groupby('accession').apply(lambda x: x.set_index('accession').to_dict('list')).to_dict()


    # create a diction with accession as key and a list of snp positions as
    # values
    snps_dict = dict()
    for index, row in df1.iterrows():
        accession = row['accession']
        positions = row['positions']

        if accession not in snps_dict:
            snps_dict[accession] = [positions]
        else:
            snps_dict[accession].append(positions)

    # generate a dict for the alignment
    f_dict = dict()
    fasta_dict = dict(fasta_iterator(alignment))
    for k, v in fasta_dict.items():
        f_dict[k.rsplit("|")[0]] = [str(base) for base in v]

    # extract flanking sequences
    flanking_seqs_dic = dict()
    for k, v in ord_dict.items():
        for index, pos in enumerate(ord_dict.get(k)['positions']):
            mutlist = [*ord_dict.get(k)['SNP'][index]]
            seqlen = len(f_dict.get(k))
            
            if 'N' in mutlist or pos >= seqlen+1:
                continue
            else:
                if ord_dict.get(k)['SNP'][index] == mutation:
                    sequence_context = ''.join(f_dict.get(k)[pos-6:pos+5]).upper()
                    print(k, pos, mutlist, mutation, f_dict.get(k)[pos], ord_dict.get(k)['SNP'][index], f_dict.get(k)[pos-6:pos], f_dict.get(k)[pos:pos+5], sequence_context, len([*sequence_context]), len(sequence_context))
                    if k not in flanking_seqs_dic:
                        flanking_seqs_dic[k] = [sequence_context]
                    else:
                        flanking_seqs_dic[k].append(sequence_context)


    # for k, v in flanking_seqs_dic.items():
    #     print(k, v)
    lineages_flanking_dic = dict()

    cols = list(columns)
    groups = df.groupby(cols)
    for name, group in groups:
        if len(cols) > 1:
            prefix = "_".join(name)
        else:
            prefix = name

        df2 = groups.get_group(name)
        d = df2.set_index("accession").T.to_dict("list")

        for k, v in d.items():
            if flanking_seqs_dic.get(k) is None:
                continue
            else:
                if prefix not in lineages_flanking_dic:
                    lineages_flanking_dic[prefix] = flanking_seqs_dic.get(k)
                else:
                    lineages_flanking_dic[prefix] += flanking_seqs_dic.get(k)

    sequences_dict = dict()
    for k, v in lineages_flanking_dic.items():
        sequences_dict[k] = sorted(list(v))
        outfile = os.path.join(outdir, str(k)) + "." + mutation + '.sequence-context.txt'
        logging.info(
            "Writing {} non-redundant sequence contexts: to {}".format(len(v), os.path.basename(outfile)))
        with open(outfile, 'w') as f_obj:
            f_obj.write(k + "," + ','.join(sorted(list(v))) + "\n")


def main():
    parser = parse_args()
    args = parser.parse_args()

    if args.alignment is not None and isinstance(args.alignment, str):
        logging.info("input alignment = {:^10}".format(os.path.basename(args.alignment), str))

    if args.snps is not None and isinstance(args.snps, str):
        logging.info("input snps metadata = {:^10}".format(os.path.basename(args.snps), str))

    outdir = os.path.abspath(args.outdir)
    mkdir(outdir)

    extract_sequences(alignment=args.alignment,
                      snps=args.snps,
                      mutation=args.mutation,
                      columns=args.columns,
                      outdir=outdir)


if __name__ == '__main__':
    main()
