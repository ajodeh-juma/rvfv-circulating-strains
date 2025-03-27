#!/usr/bin/env python3
"""
Extract flanking sequences for APOBEC sequence contexts
"""

import os
import re
import csv
import sys
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
    required_group.add_argument('--positive-strand', type=str,
                                dest="positive_strand", metavar="<str>",
                                help="path to the positive/plus/forward strand alignment file in fasta format"
                                )
    required_group.add_argument('--negative-strand', type=str,
                                dest="negative_strand", metavar="<str>",
                                help="path to the negative/minus/reverse complementary strand alignment file in fasta format"
                                )                          
    required_group.add_argument('--positive-strand-snps', type=str, 
                                dest="positive_strand_snps", metavar="<str>", 
                                help="a comma-separated file containing SNPs in positive strand, must have the columns labelled 'snps', "
                                     "'positions', 'accession'")
    required_group.add_argument('--negative-strand-snps', type=str, 
                                dest="negative_strand_snps", metavar="<str>", 
                                help="a comma-separated file containing SNPs in negative strand, must have the columns labelled 'snps', "
                                     "'positions', 'accession'")
    required_group.add_argument('--num', type=int, dest="num", metavar="<int>", default=5,
                                help="number of flanking nucleotides to extract upstream and downstream of the SNP")
    required_group.add_argument('--columns', nargs='+',
                                dest="columns", metavar="<str>",
                                choices=['year', 'host', 'lineage', 'strain_type'],
                                help="space-separated list of the type of header details to include (i.e "
                                     "year, host, lineage, strain_type")
    required_group.add_argument('--prefix', required=True, metavar='<str>',
                                dest="prefix",
                                help="prefix of the output file(s)"
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


def get_sequence_context_records(alignment, snps, num, mutation, prefix, outdir):
    """
    """

    df = pd.read_csv(snps)
    df['SNP'] = df['snps'].apply(lambda x: ''.join(re.findall(r"([A-Z]+)", str(x))))
    df1 = df[(df['SNP'] == mutation)]
    df1 = df[(df['strain_type'] == 'non-vaccine')]
    
     # create a dictionary with accession as key and list of positions and snps
    cols = ['accession', 'positions', 'SNP']
    ord_dict = df1[cols].groupby('accession').apply(lambda x: x.set_index('accession').to_dict('list')).to_dict()


    # create a dictionary with accession as key and a list of snp positions as
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


    # extract flanking sequences and associated records
    flanking_seqs_dic = dict()
    for k, v in ord_dict.items():
        for index, pos in enumerate(ord_dict.get(k)['positions']):
            pos = int(pos)
            mut = ord_dict.get(k)['SNP'][index]
            mutlist = [*mut]
            seqlen = len(f_dict.get(k))
            if 'N' in mutlist or pos >= seqlen+1:
                continue
            else:
                if (pos == snps_dict.get(k)[index] and mut == mutation):
                    snp = f_dict.get(k)[pos-1]
                    sequence_context = ''.join(f_dict.get(k)[pos-num-1:pos+num]).upper()
                    sequence_context_len = len(sequence_context)
                    if sequence_context_len == 0 or ((seqlen-pos+1) == num):
                        continue
                    else:
                        # get position of snp on the complementary sequence
                        pos_complement = (seqlen - pos) + 1

                        # get upstream nucleotides
                        upstream_contexts = [f_dict.get(k)[pos-1-num:pos-(i+1)] for i in range(num)]
                        upstream_bases = list(reversed([i[-1] for i in upstream_contexts]))
                        minus_1_position = upstream_bases[-1]
                        minus_2_position = upstream_bases[-2]
                        # minus_3_position = upstream_bases[-3]


                        # get downstream nucleotides
                        downstream_contexts = [f_dict.get(k)[pos-1+(num-i):pos+num] for i in range(num)]
                        downstream_bases = list(reversed([i[0] for i in downstream_contexts]))
                        plus_1_position = downstream_bases[0]
                        plus_2_position = downstream_bases[1]
                        # plus_3_position = downstream_bases[2]

                        dinucleotide = ''.join([minus_2_position, minus_1_position])
                        hotspot_di = ''.join([minus_1_position, snp])
                        hotspot_tri = ''.join([minus_2_position, minus_1_position, snp])

                        # print(k, pos, mut, mutation, seqlen, snp, sequence_context, pos_complement, minus_1_position, minus_2_position, minus_3_position, plus_1_position, plus_2_position, plus_3_position, dinucleotide, hotspot_di, hotspot_tri)

                        flanking_seqs_dic[k + ' ' + str(pos)] = [k + ' ' + str(pos), k, pos, pos_complement, snp, mut, sequence_context, minus_1_position, minus_2_position, plus_1_position, plus_2_position, dinucleotide, hotspot_di, hotspot_tri, seqlen]

    col_names=['accession_pos', 'accession', 'position', 'position_complement', 'snp', 'mutation', 'sequence_context', 'minus_1_position', 'minus_2_position', 'plus_1_position', 'plus_2_position', 'dinucleotide', 'hotspot_di', 'hotspot_tri', 'length']
    sequence_contexts_df = pd.DataFrame.from_dict(flanking_seqs_dic, orient='index', columns=col_names)
    sequence_contexts_df.reset_index(drop=True, inplace=True)

    outfile = os.path.join(outdir, prefix + '-' + mutation + '-' + 'sequence-contexts-with-metadata.csv')
    sequence_contexts_df.to_csv(outfile, index=False)
    return sequence_contexts_df, df


def identify_apobecs_on_strands(plus_df, minus_df, snps_df, columns, prefix, outdir):
    """
    """


    outfile = os.path.join(outdir, prefix + '-' + 'sequence-contexts-with-metadata.csv')
    

    # merge dataframes to get potential apobec mutated sites on plus and minus strands
    merged_df = pd.merge(plus_df, minus_df, left_on=['accession', 'position', 'length'], right_on=['accession', 'position_complement', 'length'], how="inner")

    # renames columns
    new_colnames = {'accession_pos_x': 'accession_plus_strand', 
                    'position_x': 'position_plus_strand', 
                    'position_complement_x': 'position_complement_minus_strand',
                    'snp_x': 'snp_plus_strand', 
                    'mutation_x': 'mutation_plus_strand', 
                    'sequence_context_x': 'sequence_context_plus_strand', 
                    'minus_1_position_x': 'minus_1_pos_plus_strand',
                    'minus_2_position_x': 'minus_2_pos_plus_strand', 
                    'plus_1_position_x': 'plus_1_pos_plus_strand',
                    'plus_2_position_x': 'plus_2_pos_plus_strand', 
                    'dinucleotide_x': 'dinucleotide_plus_strand',
                    'hotspot_di_x': 'hotspot_di_plus_strand',
                    'hotspot_tri_x': 'hotspot_tri_plus_strand',
                    'accession_pos_y': 'accession_minus_strand', 
                    'position_y': 'position_minus_strand', 
                    'position_complement_y': 'position_complement_plus_strand', 
                    'snp_y': 'snp_minus_strand',
                    'mutation_y': 'mutation_minus_strand', 
                    'sequence_context_y': 'sequence_context_minus_strand',
                    'minus_1_position_y': 'minus_1_pos_minus_strand',
                    'minus_2_position_y': 'minus_2_pos_minus_strand', 
                    'plus_1_position_y': 'plus_1_pos_minus_strand',
                    'plus_2_position_y': 'plus_2_pos_minus_strand', 
                    'dinucleotide_y': 'dinucleotide_minus_strand',
                    'hotspot_di_y': 'hotspot_di_minus_strand',
                    'hotspot_tri_y': 'hotspot_tri_minus_strand',}

    
    merged_df = merged_df.rename(columns=new_colnames)
    merged_df = merged_df[['accession', 
                            'position_plus_strand', 
                            'snp_plus_strand', 
                            'mutation_plus_strand', 
                            'sequence_context_plus_strand', 
                            'minus_1_pos_plus_strand', 
                            'minus_2_pos_plus_strand', 
                            'plus_1_pos_plus_strand', 
                            'plus_2_pos_plus_strand', 
                            'dinucleotide_plus_strand',
                            'hotspot_di_plus_strand',
                            'hotspot_tri_plus_strand',
                            'position_minus_strand', 
                            'snp_minus_strand',
                            'mutation_minus_strand', 
                            'sequence_context_minus_strand',
                            'minus_1_pos_minus_strand',
                            'minus_2_pos_minus_strand',
                            'plus_1_pos_minus_strand',
                            'plus_2_pos_minus_strand',
                            'dinucleotide_minus_strand',
                            'hotspot_di_minus_strand',
                            'hotspot_tri_minus_strand',
                            'length']]

    snps_df = snps_df[['accession', 'taxa', 'host', 'lineage', 'strain_type', 'year']].drop_duplicates()
    print(len(snps_df))

    host_groups = snps_df.groupby(['host'])
    for host, group in host_groups:
        df1_host = host_groups.get_group(host)
        print(host, len(df1_host))


    df3 = pd.merge(merged_df, snps_df, left_on=['accession'], right_on=['accession'], how="left")
    logging.info(
            "Writing {} all sequence contexts: to {}".format(len(df3), os.path.basename(outfile)))
    df3.to_csv(outfile, index=False)
    
    
    plus_dict = dict()
    minus_dict = dict()
    for index, row in df3.iterrows():
        k = row['accession']

        if k not in plus_dict:
            plus_dict[k] = [row['sequence_context_plus_strand']]
        else:
            plus_dict[k].append(row['sequence_context_plus_strand'])


        if k not in minus_dict:
            minus_dict[k] = [row['sequence_context_minus_strand']]
        else:
            minus_dict[k].append(row['sequence_context_minus_strand'])

    return plus_dict, minus_dict, df3


def write_records(dic, snps_df, columns, mutation, prefix, outdir):
    """
    """

    snps_df = snps_df[['accession', 'taxa', 'host', 'lineage', 'strain_type', 'year']].drop_duplicates()

    by_columns_dict = dict()

    cols = list(columns)
    groups = snps_df.groupby(cols)
    for name, group in groups:
        if len(cols) > 1:
            selected_cols = "_".join(name)
        else:
            selected_cols = name

        df2 = groups.get_group(name)
        d = df2.set_index("accession").T.to_dict("list")

        for k, v in d.items():
            if dic.get(k) is None:
                continue
            else:
                if selected_cols not in by_columns_dict:
                    by_columns_dict[selected_cols] = dic.get(k)
                else:
                    by_columns_dict[selected_cols] += dic.get(k)

    sequences_dict = dict()
    for k, v in by_columns_dict.items():
        # sequences_dict[k] = sorted(list(v))
        outfile = os.path.join(outdir, str(k)) + "." + mutation + "." + prefix + ".sequence-context.txt"
        logging.info(
            "Writing {} non-redundant sequence contexts: to {}".format(len(v), os.path.basename(outfile)))
        with open(outfile, 'w') as f_obj:
            f_obj.write(k + "," + ','.join(list(v)) + "\n")



def classify_apobecs(df, prefix, outdir):
    """
    """

    # get number of sequences per host
    groups = df.groupby(['host'])
    for name, group in groups:
        df2 = groups.get_group(name)
        df2_sel = df2[['accession', 'lineage', 'host', 'position_plus_strand', 'position_minus_strand', 'sequence_context_plus_strand', 'sequence_context_minus_strand']]
        # print(df2_sel.sort_values(by=['position_plus_strand'], ascending=True   ))

    # apobec3g
    df1 = df[(df.dinucleotide_plus_strand == 'gg') & (df.dinucleotide_minus_strand == 'cc')]
    df1['apobec_type'] = 'A3G'
    filename = 'A3G' + '.' + prefix + ".edited.sites.txt"
    outfile = os.path.join(outdir, filename)
    df1.to_csv(outfile, index=False)
    write_classified_apobecs(df=df1, prefix=prefix, outdir=outdir)


    # apobec3a/apobec3b
    df2 = df[(df.minus_1_pos_plus_strand == 'a') & (df.snp_plus_strand == 'a') & (df.dinucleotide_plus_strand != 'aa') & (df.dinucleotide_plus_strand != 'gg')
            & (df.minus_1_pos_minus_strand == 't') & (df.snp_minus_strand == 't') & (df.dinucleotide_minus_strand != 'cc') & (df.dinucleotide_minus_strand != 'aa')]
    df2['apobec_type'] = 'A3A_A3B'
    filename = 'A3A_A3B' + '.' + prefix + ".edited.sites.txt"
    outfile = os.path.join(outdir, filename)
    df2.to_csv(outfile, index=False)
    write_classified_apobecs(df=df2, prefix=prefix, outdir=outdir)



    # apobec3c/apobec3f
    df3 = df[(df.dinucleotide_plus_strand == 'aa') & (df.snp_plus_strand == 'a') & (df.dinucleotide_minus_strand == 'tt') & (df.snp_minus_strand == 't')]
    df3['apobec_type'] = 'A3C_A3F'
    filename = 'A3C_A3F' + '.' + prefix + ".edited.sites.txt"
    outfile = os.path.join(outdir, filename)
    df3.to_csv(outfile, index=False)
    write_classified_apobecs(df=df3, prefix=prefix, outdir=outdir)


def write_classified_apobecs(df, prefix, outdir):
    """
    """

    d_plus = dict()
    d_minus = dict()
    if len(df) > 0:
        for index, row in df.iterrows():
            # print(index, row['apobec_type'], row['sequence_context_plus_strand'], row['sequence_context_minus_strand'])
            if not row['apobec_type'] in d_plus:
                d_plus[row['apobec_type']] = [row['sequence_context_plus_strand']]
                d_minus[row['apobec_type']] = [row['sequence_context_minus_strand']]
            else:
                d_plus[row['apobec_type']].append(row['sequence_context_plus_strand'])
                d_minus[row['apobec_type']].append(row['sequence_context_minus_strand'])

    for k, v in d_plus.items():
        filename = str(k) + '.' + "GA" + '.' + prefix + ".context.txt"
        outfile = os.path.join(outdir, filename)
        with open(outfile, 'w') as f_obj:
            f_obj.write(k + "," + ','.join(v) + "\n")

    for k, v in d_minus.items():
        filename = str(k) + '.' + "CT" + '.' + prefix + ".context.txt"
        outfile = os.path.join(outdir, filename)
        with open(outfile, 'w') as f_obj:
            f_obj.write(k + "," + ','.join(v) + "\n")



def main():
    parser = parse_args()
    args = parser.parse_args()

    if args.positive_strand is not None and isinstance(args.positive_strand, str):
        logging.info("input positive strand alignment = {:^10}".format(args.positive_strand, str))

    if args.positive_strand_snps is not None and isinstance(args.positive_strand_snps, str):
        logging.info("input positive snps = {:^10}".format(args.positive_strand_snps, str))

    if args.negative_strand is not None and isinstance(args.negative_strand, str):
        logging.info("input negative strand alignment = {:^10}".format(args.negative_strand, str))

    if args.negative_strand_snps is not None and isinstance(args.negative_strand_snps, str):
        logging.info("input negative snps = {:^10}".format(args.negative_strand_snps, str))


    if args.num < 2:
        logging.info("number must be greater or equal to 2")
        sys.exit(1)

    outdir = os.path.abspath(args.outdir)
    mkdir(outdir)


    plusdf, plus_snps_df = get_sequence_context_records(alignment=args.positive_strand, 
                    snps=args.positive_strand_snps,
                    num=args.num,
                    mutation='GA',
                    prefix=args.prefix,
                    outdir=outdir)

    minusdf, minus_snps_df = get_sequence_context_records(alignment=args.negative_strand, 
                    snps=args.negative_strand_snps,
                    num=args.num,
                    mutation='CT',
                    prefix=args.prefix,
                    outdir=outdir)

    plus_dict, minus_dict, merged = identify_apobecs_on_strands(plus_df=plusdf, 
                                                        minus_df=minusdf, 
                                                        snps_df=minus_snps_df, 
                                                        columns=args.columns,
                                                        prefix=args.prefix,
                                                        outdir=outdir)

    write_records(dic=plus_dict, 
                    snps_df=plus_snps_df, 
                    columns=args.columns, 
                    mutation='GA',
                    prefix=args.prefix,
                    outdir=outdir)

    write_records(dic=minus_dict, 
                    snps_df=minus_snps_df, 
                    columns=args.columns, 
                    mutation='CT',
                    prefix=args.prefix,
                    outdir=outdir)

    classify_apobecs(merged, prefix=args.prefix, outdir=outdir)


if __name__ == '__main__':
    main()
