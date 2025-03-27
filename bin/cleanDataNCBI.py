#!/usr/bin/env python3

"""
clean data downloaded from NCBI using correct sequence information from literature or other reports and filter out
potential reassortants and other outliers
"""

import os
import logging
import argparse
import pycountry
import textwrap
from textwrap import dedent
import pandas as pd
from utils import fasta_iterator

log_level = logging.DEBUG
logging.basicConfig(level=log_level, format='[%(asctime)s] - %(''levelname)s - [%(funcName)s] - %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S %p')


def parse_args():
    # arguments parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        argument_default=argparse.SUPPRESS, prefix_chars='--',
        description=dedent('''clean data downloaded from NCBI using correct sequence information from literature or other reports and filter out 
potential reassortants and other outliers''')
    )
    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument('--fasta', type=str, dest="fasta", metavar="<FILE>",
                                help="FASTA file format having sequences obtained from NCBI with the header being the "
                                     "accession number "
                                )
    required_group.add_argument('--metadata', type=str, dest="metadata", metavar="<FILE>",
                                help="comma-separated values (csv) file having the sequence metadata and columns "
                                     "'accession, organism, length, strain, country, location, host and date'"
                                )
    required_group.add_argument('--strains-metadata', type=str, dest="strains_metadata", metavar='<FILE>',
                                help="comma-separated values (csv) file having the virus strains metadata and columns "
                                     "'strain, country, location, host and date'. The strain column should match "
                                     "strain values in the metadata file "
                                )
    required_group.add_argument('--prefix', required=True, metavar='<str>', dest="prefix",
                                help="prefix of the output file(s)"
                                )
    required_group.add_argument('--min-length', required=True, type=int,
                                dest="min_length", metavar="<INT>",
                                help="Minimum length of the sequence record"
                                )
    parser.add_argument('--outDir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default="./")
    return parser



host_dict = {
        'Aedes circumluteolus': 'Mosquito',
        'Aedes circumluteolus (mosquito)': 'Mosquito',
        'Aedes cuminsi': 'Mosquito',
        'Aedes cumminsii (mosquito)': 'Mosquito',
        'Aedes dalzieli (mosquito)': 'Mosquito',
        'Aedes macintoshi': 'Mosquito',
        'Aedes mcintoshi (mosquito)': 'Mosquito',
        'Aedes ochraceus': 'Mosquito',
        'Aedes palpalis': 'Mosquito',
        'Aedes sp. (mosquito)': 'Mosquito',
        'Aedes taeniorhynchus': 'Mosquito',
        'Aedes vexans arabiensis': 'Mosquito',
        'African buffalo': 'Buffalo',
        'Bos taurus': 'Cow',
        'Buffalo': 'Buffalo',
        'Cephalophorus callipygus': 'Antelope',
        'Culex pipiens': 'Mosquito',
        'Culex zombaensis': 'Mosquito',
        'Culex zombaensis (mosquito)': 'Mosquito',
        'Eretmapodites quinquevittatus': 'Mosquito',
        'Eretmapodites quinquevittatus (mosquito)': 'Mosquito',
        'Hipposideros caffer (bat)': 'Bat',
        'Homo sapiens': 'Human',
        'Mansonia africana': 'Mosquito',
        'Micropteropus pusillus (bat)': 'Bat',
        'Micropterus pusillus (bat)': 'Bat',
        'Philantomba monticola': 'Antelope',
        'Syncerus caffer (buffalo)': 'Buffalo',
        'Vero cells': 'Human',
        'animal': 'Sheep',
        'bovine': 'Cow',
        'buffalo': 'Buffalo',
        'caprine': 'Goat',
        'cow': 'Cow',
        'mosquito': 'Mosquito',
        'ovine': 'Sheep',
        'plaque pick of RVFV 74HB59': 'Human',
        'sheep': 'Sheep',
        'small plaque pick of ZH-501': 'Human',
        'springbok': 'Antelope',
        'Aedes macintoshi': 'Mosquito',
        'Aedes ochraceus': 'Mosquito',
        'Aedes palpalis': 'Mosquito',
        'Aedes taeniorhyncus': 'Mosquito',
        'Aedes vexans arabiensis': 'Mosquito',
        'African buffalo': 'Buffalo',
        'Bovine calf': 'Cow',
        'Buffalo calf': 'Cow',
        'Culex pipiens' : 'Mosquito',
        'Homo sapiens': 'Human',
        'Human': 'Human',
        'Human sera': 'Human',
        'Mosquito': 'Mosquito',
        'Ovine': 'Sheep',
        'bat': 'Bat',
        'goat': 'Goat',
        'homo sapiens': 'Human',
        'Bovine': 'Cow',
        'human': 'Human',
        'human sera': 'Human',
    }

def read_fasta(fasta):
    """
    read fasta file and convert to dictionary
    :param fasta:
    :return:
    """

    seq_dict = dict((x[0].split()[0], x[1]) for x in fasta_iterator(fasta_name=fasta))
    return seq_dict


def read_meta(metadata_fn):
    """
    read metadata csv file
    :param metadata_fn:
    :return:
    """

    df = pd.read_csv(metadata_fn, sep=",")
    return df


def rename_host(host):
    """

    @param host:
    @return:
    """

    if 'sapiens' in host or 'traveller' in host or 'patient' in host or 'serum' in host or 'Human' in host or 'cerebrospinal fluid' in host:
        return 'Human'
    elif 'Aedes' in host or 'aedes' in host or 'Anopheles' in host or 'Culex' in host:
        return 'Mosquito'
    elif 'mouse' in host:
        return 'Mouse'
    elif 'Macaca' in host:
        return 'Monkey'
    elif 'micropterus' or 'hipposideros' or 'eretmapodites' or 'micropteropus' in host:
        return 'Bat'
    else:
        return host


def clean_data(seq_dict, df1, df2, min_length, prefix, outdir):
    """

    :param seq_dict:
    :param df1:
    :param df2:
    :param min_length:
    :param prefix:
    :param outdir:
    :return:
    """

    # filter on length
    df1_on_len = df1[(df1.length.astype(int) > min_length)]
    to_exclude = len(df1) - len(df1_on_len)
    logging.info("{} sequences with sequence lengths < {} excluded".format(to_exclude, min_length))

    # add correct metadata information
    df3 = pd.merge(df1, df2, how="inner", on="strain")
    df4 = df3[['accession', 'organism', 'length',
               'strain', 'country_y', 'location_y', 'host_y', 'date_y']].rename(columns={'country_y': 'country',
                                                                                         'location_y': 'location',
                                                                                         'host_y': 'host',
                                                                                         'date_y': 'date'})
    print(df4)
    print(list(set(df4.host)))
    # fill NAs with empty strings
    df4 = df4.fillna('')

    # renames hosts
    df4.host = df4.host.map(host_dict)
    df4.host = df4.host.str.lower()
    print(list(set(df4.host)))

    # get unshared accessions
    non_common_accessions = set(list(df4.accession)) ^ set(list(df1_on_len.accession))
    print(non_common_accessions)
    logging.info("{} sequences without records".format(len(non_common_accessions)))
    # print(df1_on_len[df1_on_len.accession.isin(non_common_accessions)])
    logging.info("{} clean sequences found".format(len(df4)))

    # write metadata
    outfile_1 = os.path.join(outdir, prefix) + '.csv'
    df4.to_csv(outfile_1, sep=",", index=False)

    # create a dictionary from clean metadata dataframe
    d = df4.set_index('accession').T.to_dict('list')

    # write fasta

    # use module pycountry to create a dictionary of counties and their 3-letter codes/symbols
    countries = dict()
    t = list(pycountry.countries)
    for country in t:
        countries[country.name.split(',')[0]] = country.alpha_3

    outfile_2 = os.path.join(outdir, prefix) + '.fasta'
    with open(outfile_2, 'w') as f_out:
        for accession, record in d.items():
            # print(accession, record)
            # new header format >accession|host|country|location|date
            header = accession + '|' + record[5] + '|' + countries.get(record[3]) + '|' + record[4] + '|' + record[6]
            sequence = seq_dict.get(accession.rstrip('>'))
            f_out.write((">{}\n{}\n".format(header, textwrap.fill(sequence, width=80))))

    return outfile_1, outfile_2


def main():
    parser = parse_args()
    args = parser.parse_args()

    d = read_fasta(fasta=args.fasta)
    df1 = read_meta(metadata_fn=args.metadata)
    df2 = read_meta(metadata_fn=args.strains_metadata)

    # create outdir if not exists
    outdir = os.path.abspath(args.outdir)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    clean_data(seq_dict=d, df1=df1, df2=df2, min_length=args.min_length, prefix=args.prefix, outdir=args.outdir)


if __name__ == '__main__':
    main()
