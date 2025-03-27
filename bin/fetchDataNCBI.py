#!/usr/bin/env python3

"""
fetch viral sequence data and related metadata from ncbi
"""

import os
import re
import sys
import logging
import argparse
import textwrap
from textwrap import dedent
import pandas as pd
from collections import OrderedDict

import pycountry
from Bio import SeqIO
from Bio import Entrez

from utils import fasta_iterator
from utils import str2bool

log_level = logging.DEBUG
logging.basicConfig(level=log_level, format='[%(asctime)s] - %(''levelname)s - [%(funcName)s] - %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S %p')


def parse_args():
    # arguments parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        argument_default=argparse.SUPPRESS, prefix_chars='--',
        description=dedent(
            '''fetch viral sequence data and metadata from ncbi given a text file with accession numbers''')
    )
    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument('--accessions', type=str, dest="accessions", metavar="<FILE>",
                                help="comma-separated values (csv) file having a column with accessions labelled "
                                     "as 'accession'"
                                )
    required_group.add_argument('--database', type=str, metavar='<str>', default="nucleotide",
                                choices=["nucleotide", "gene", "protein"], help="NCBI database name"
                                )
    parser.add_argument('--outDir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default="./"
                        )
    parser.add_argument('--batch-size', type=int, required=False, metavar='batch_size', default=100,
                        help='number of accessions to process per request'
                        )
    parser.add_argument('--email', type=str, metavar="<email>", dest="email", required=False,
                        default="someone@gmail.com", help="email address to access Entrez"
                        )
    parser.add_argument('--cleanup', metavar="<boolean>", type=str2bool, default=False,
                        dest='cleanup',
                        help=dedent('''delete the records (genbank (.gbk) and fasta) output files once downloaded''')
                        )
    return parser


def read_accessions(accessions):
    """

    :param accessions:
    :return:
    """

    # read the csv file
    lineages_df = pd.read_csv(accessions)

    # filter out rows with no accessions
    lineages_df = lineages_df[lineages_df['accession'].notnull()]

    # drop duplicates
    lineages_df = lineages_df.drop_duplicates(subset=['accession'], keep='first', inplace=False)
    accessions = list(lineages_df.accession)
    return accessions


def accessions_to_gb(accessions, email, database, batch_size, ret_max=10 ** 9):
    """

    :param accessions:
    :param database:
    :param email:
    :param batch_size:
    :param ret_max:
    :return:
    """

    Entrez.email = repr(email)  # email address

    def batch(sequence, size):
        """

        :param sequence:
        :param size:
        :return:
        """
        length = len(accessions)
        for start in range(0, length, size):
            yield sequence[start:min(start + size, length)]

    def extract_records(records_handle):
        """

        :param records_handle:
        :return:
        """
        logging.info("getting records for accession")
        buffer = []
        for line in records_handle:
            if line.startswith("LOCUS") and buffer:
                # yield accession number and record
                yield buffer[0].split()[1], "".join(buffer)
                buffer = [line]
            else:
                buffer.append(line)
        yield buffer[0].split()[1], "".join(buffer)

    def process_batch(accessions_batch):
        """

        :param accessions_batch:
        :return:
        """
        # get GI for query accessions
        logging.info("getting GI accessions")
        query = " ".join(accessions_batch)
        query_handle = Entrez.esearch(db=database, term=query, retmax=ret_max)
        gi_list = Entrez.read(query_handle, validate=False)['IdList']

        # get GB files
        logging.info("getting the genbank file format")
        search_handle = Entrez.epost(db=database, id=",".join(gi_list))
        search_results = Entrez.read(search_handle, validate=False)
        webenv, query_key = search_results["WebEnv"], search_results["QueryKey"]
        records_handle = Entrez.efetch(db=database, rettype="gb", retmax=batch_size,
                                       webenv=webenv, query_key=query_key)

        yield from extract_records(records_handle)

    accession_batches = batch(accessions, batch_size)
    # print(accession_batches)
    for acc_batch in accession_batches:
        print(acc_batch)
        yield from process_batch(acc_batch)


def write_genbank(accession, record, outdir):
    """
    :param accession:
    :param record:
    :param outdir:
    :return:
    """

    # write genbank record
    gbk = os.path.join(outdir, accession + '.gbk')
    with open(gbk, "w") as f_obj:
        print(record, file=f_obj)
    return gbk


def write_fasta(accession, genbank, outdir):
    """

    :param accession:
    :param genbank:
    :param outdir:
    :return:
    """
    # write fasta
    fasta = os.path.join(outdir, accession + '.fasta')
    records = SeqIO.parse(genbank, "genbank")
    with open(fasta, 'w') as f_obj:
        for record in records:
            f_obj.write((">{}\n{}\n".format(accession, textwrap.fill(str(record.seq), width=80))))
    return fasta


def write_metadata(accession, dataframe, outdir):
    """

    :param accession:
    :param genbank:
    :param outdir:
    :return:
    """
    # write metadata
    metadata_csv = os.path.join(outdir, accession + '.csv')
    df = dataframe[['accession', 'organism', 'length', 'strain', 'country', 'location', 'host', 'date']]
    df.to_csv(metadata_csv, sep=",", index=False)
    return metadata_csv


def write_output(metadata_fns, fasta_fns, outdir):
    """

    :param metadata_fns:
    :param fasta_fns:
    :param outdir:
    @return:
    """

    fasta = os.path.join(outdir, 'sequences.fasta')
    metadata = os.path.join(outdir, 'sequences.csv')

    # write all fasta sequences into a single file
    with open(fasta, 'wt') as fh:
        for fn in fasta_fns:
            if os.path.isfile(fn):
                seq_dict = dict((x[0].split()[0], x[1]) for x in fasta_iterator(fasta_name=fn))
                for seq_id, seq in seq_dict.items():
                    # logging.info("writing {} sequences into {}".format(fn, fasta))
                    fh.write((">{}\n{}\n".format(seq_id, textwrap.fill(seq, width=80))))
    logging.info("wrote {} FASTA sequences to {}".format(len(fasta_fns), fasta))
    # write all metadata into a single file
    with open(metadata, 'wt') as f_out:
        f_out.write("accession,organism,length,strain,country,location,host,date")
        f_out.write("\n")
        for fn in metadata_fns:
            if os.path.isfile(fn):
                for line in open(fn):
                    if line.startswith('accession'):
                        continue
                    else:
                        f_out.write("\t".join(line.strip().split('\t')))
                        f_out.write("\n")
    logging.info("wrote {} metadata records to {}".format(len(metadata_fns), metadata))
    return metadata, fasta


class GenBankParser:

    def __init__(self, gbk):
        """

        @param gbk:
        """
        self.records = SeqIO.parse(gbk, "genbank")  # parse records from the genbank file
        self.organism = ""
        self.accession = ""
        self.sequence = ""
        self.sequence_length = ""
        self.annotations = ""
        self.strain = ""
        self.country = ""
        self.location = ""
        self.host = ""
        self.date = ""
        self.pubmed_id = ""

    def extract_records(self):
        """
        @return:
        """

        # use module pycountry to create a dictionary of counties and their 3-letter codes/symbols
        # countries = dict()
        # t = list(pycountry.countries)
        # for country in t:
        #     countries[country.name.split(',')[0]] = country.alpha_3

        # create a dictionary with accession as key and records as values
        d = OrderedDict()
        for record in self.records:
            source = record.features[0]
            self.accession = record.id.rsplit(".")[0]
            self.organism = record.annotations['organism']
            self.sequence_length = len(record.seq)
            self.sequence = str(record.seq)
            self.pubmed_id = record.annotations['references'][0].pubmed_id

            # get strain
            try:
                self.strain = ''.join(source.qualifiers['strain'])
            except KeyError:
                self.strain = ''
                if KeyError:
                    try:
                        self.strain = ''.join(source.qualifiers['isolate'])
                    except KeyError:
                        self.strain = ''
                        if KeyError:
                            try:
                                self.strain = ''.join(source.qualifiers['clone'])
                            except KeyError:
                                self.strain = ''

            # get country
            try:
                self.country = ''.join(source.qualifiers['geo_loc_name']).rsplit(':')
                #print(self.country, len(self.country))
                if len(self.country) == 1:
                    self.country = self.country[0]
                    # print(self.country)
                if len(self.country) == 2:
                    self.country, self.location = self.country[0], self.country[1].lstrip()
                    if len(self.location.split()) == 3 and self.location.split()[2] == 'State':
                        self.location = ' '.join(self.location.split()[0:2])
                    elif len(self.location.split()) == 3 and self.location.split()[2] == 'district':
                        self.location = ' '.join(self.location.split()[0:2])
                    elif len(self.location.split()) == 2:
                        self.location = self.location.split()[0]
                    else:
                        self.location = self.location.split()[0].replace(",", "")

            except KeyError:
                self.country = ''
                self.location = ''
                # if KeyError:
                #     try:
                        
                #     except KeyError:
                #         self.location = ''
            # get host/isolation source
            try:
                self.host = ''.join(source.qualifiers['host'])
            except KeyError:
                self.host = ''
                if KeyError:
                    try:
                        self.host = ''.join(source.qualifiers['isolation_source'])
                    except KeyError:
                        self.host = ''

            # get collection date
            try:
                self.date = ''.join(source.qualifiers['collection_date'])
            except KeyError:
                self.date = ''
                if KeyError:
                    try:
                        self.date = re.findall(r"[0-9]{4,11}", self.date)
                        self.date = ''.join(self.date)
                    except KeyError:
                        self.date = ''

            d[self.accession] = [self.accession, self.organism, self.sequence_length, self.strain, self.country,
                                 self.location, self.host, self.date, self.pubmed_id, self.sequence]

        columns = ['accession', 'organism', 'length', 'strain', 'country',
                   'location', 'host', 'date', 'pubmed_id', 'sequence']

        cols_d = dict()
        for i, col in enumerate(columns):
            cols_d[i] = col

        df = pd.DataFrame.from_dict(d, orient='index')
        df = df.rename(columns=cols_d).reindex()

        # change date column to datetime datatype
        df['date'] = pd.to_datetime(df['date'])
        return df


def main():
    parser = parse_args()
    args = parser.parse_args()

    # check the arguments
    if not os.path.isfile(args.accessions):
        logging.error("File {} does not exist")
        sys.exit(1)
    else:
        accessions_df = read_accessions(args.accessions)

    # create outdir if not exists
    outdir = os.path.abspath(args.outdir)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # enumerate Genbank records
    fasta_files = []
    metadata_files = []
    gbk_files = []
    for accession, record in accessions_to_gb(accessions=accessions_df,
                                              email=args.email,
                                              database=args.database,
                                              batch_size=args.batch_size,
                                              ret_max=args.batch_size):
        gbk_fn = write_genbank(accession=accession, record=record, outdir=outdir)
        fasta_fn = write_fasta(accession=accession, genbank=gbk_fn, outdir=outdir)
        GBObject = GenBankParser(gbk=gbk_fn)
        dataframe = GBObject.extract_records()
        metadata_fn = write_metadata(accession=accession, dataframe=dataframe, outdir=outdir)

        metadata_files.append(metadata_fn)
        fasta_files.append(fasta_fn)
        gbk_files.append(gbk_fn)

    write_output(fasta_fns=fasta_files, metadata_fns=metadata_files, outdir=args.outdir)

    if args.cleanup:
        [os.remove(fn) for fn in fasta_files]
        [os.remove(fn) for fn in metadata_files]
        [os.remove(fn) for fn in gbk_files]


if __name__ == '__main__':
    main()
