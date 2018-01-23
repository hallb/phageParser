#!/usr/bin/env python

import argparse
import os
import pickle

import pandas
import requests
from Bio import Entrez, SeqIO
from lxml import html, etree
from tqdm import tqdm

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'phageAPI.settings')
import django

django.setup()

from util.acc import read_accession_file
from util.prunedict import prune_dict
from util import fetch
from restapi.models import (
    Organism,
    Spacer,
    Repeat,
    LocusSpacerRepeat,
    AntiCRISPR,
    Locus
)


def populate_organism(data_dir):
    def add_organism(name, accession):
        # get the object, this also checks for duplicates
        o, created = Organism.objects.get_or_create(
            name=name, accession=accession)
        return o

    def merge_acc_names(accession_list):
        acc_name_dict = {}
        db = "nuccore"
        # Doing batches of 200 to make sure requests to NCBI are not too big
        for i in range(0, len(accession_list), 200):
            j = i + 200

            result_handle = Entrez.efetch(
                db=db, rettype="gb", id=accession_list[i:j])

            # Populate result per organism name
            records = SeqIO.parse(result_handle, 'genbank')
            for record in tqdm(records):
                # Using NCBI name, which should match accession number passed
                acc_name_dict[record.name] = record.annotations['organism']
        return acc_name_dict

    with open(os.path.join(data_dir, 'bac_accession_list.txt')) as f:
        accession_list = list(read_accession_file(f))

    acc_name_dict = merge_acc_names(accession_list)
    for acc in acc_name_dict:
        add_organism(name=acc_name_dict[acc], accession=acc)


class FastaFetcher:
    def __init__(self, url, filename, limit=-1):
        self.limit = limit
        self.url = url
        self.filename = filename

    def fasta_records(self):
        fetch.fetch(self.filename, self.url)
        result = list(SeqIO.parse(self.filename, 'fasta'))
        if self.limit > -1:
            return result[:self.limit]
        return result


def initialize_spacer_and_repeat_fetchers(data_dir, limit):
    spacer_fetcher = FastaFetcher(
        'http://crispr.i2bc.paris-saclay.fr/crispr/BLAST/Spacer/Spacerdatabase',
        os.path.join(data_dir, 'spacerdatabase.txt'),
        limit=limit)
    repeat_fetcher = FastaFetcher(
        'http://crispr.i2bc.paris-saclay.fr/crispr/BLAST/DR/DRdatabase',
        os.path.join(data_dir, 'repeatdatabase.txt'),
        limit=limit)
    return spacer_fetcher, repeat_fetcher


def repeat_records_to_dict(repeat_records):
    repeat_dict = {}
    for repeat in repeat_records:
        accession_list = repeat.name.split('|')
        sequence = str(repeat.seq)
        for accession in accession_list:
            repeat_dict[accession] = {'RepeatSeq': sequence}
    return repeat_dict


def add_spacers_to_dict(gene_dict, spacer_records):
    for spacer in spacer_records:
        accession_list = spacer.name.split('|')
        sequence = str(spacer.seq)
        for accession in accession_list:
            accession_elements = accession.split('_')
            order = accession_elements[-1]
            accession_id = '_'.join(accession_elements[:-1])
            try:
                if 'Spacers' in gene_dict[accession_id]:
                    gene_dict[accession_id]['Spacers'][order] = sequence
                else:
                    gene_dict[accession_id]['Spacers'] = {order: sequence}
            except KeyError:
                print('Error on accession id:  %s' % accession_id)
    return gene_dict


def addpositionstodict(gendict):
    print("Downloading position information from web...")
    for accidwithloc in tqdm(gendict):
        if 'Start' in gendict[accidwithloc]:
            continue
        accid = '_'.join(accidwithloc.split('_')[:-1])
        url = ('http://crispr.i2bc.paris-saclay.fr/crispr/crispr_db.php?'
               'checked%5B%5D={}'.format(accid))
        page = requests.get(url)
        htmltable = html.fromstring(page.content).xpath(
            "//table[normalize-space(@class)='primary_table']")[1]
        strtable = etree.tostring(htmltable)
        # converts to pandas df and then to numpy array then drop titles
        arrtable = pandas.read_html(strtable)[0].as_matrix()[2:]
        for row in arrtable:
            if row[0] in gendict:
                gendict[row[0]]['Start'] = row[2]
                gendict[row[0]]['Stop'] = row[3]
            else:
                if row[1] != 'questionable':
                    print("Can't find %s in local files" % row[0])
    return gendict


def populate_fromlocus(locid, locdict):
    accid = '_'.join(locid.split('_')[:-1])
    organismset = Organism.objects.filter(accession=accid)
    if not organismset.exists():
        print('Organism with accid %s not found in db' % accid)
        return
    organism = organismset[0]
    repeat, _ = Repeat.objects.get_or_create(sequence=locdict['RepeatSeq'])
    loc_start = int(locdict['Start'])
    loc_end = int(locdict['Stop'])
    locus, _ = Locus.objects.get_or_create(
        organism=organism,
        genomic_start=loc_start,
        genomic_end=loc_end
    )
    spacers = locdict['Spacers']
    for order in sorted(spacers):
        spacer, _ = Spacer.objects.get_or_create(sequence=spacers[order])
        order = int(order)
        lsr, _ = LocusSpacerRepeat.objects.get_or_create(
            locus=locus,
            spacer=spacer,
            repeat=repeat,
            order=order
        )
        spacer.save()
        lsr.save()
    locus.save()
    repeat.save()
    organism.save()


def populate_lsrpair(data_dir, limit):
    print('Downloading files and gathering online data.')
    spacer_fetcher, repeat_fetcher = initialize_spacer_and_repeat_fetchers(data_dir, limit)
    gendict = prune_dict(
        addpositionstodict(
            add_spacers_to_dict(
                repeat_records_to_dict(repeat_fetcher.fasta_records()),
                spacer_fetcher.fasta_records())))
    with open(os.path.join(data_dir, 'gene_dict.pickle'), 'wb') as f:
        pickle.dump(gendict, f, protocol=pickle.HIGHEST_PROTOCOL)

    print('Created dictionary and dumped data to gene_dict.pickle')
    print('Populating Spacer, Repeat, SpacerRepeatPair, OrganismSpacerRepeatPair tables')
    for locid in tqdm(gendict):
        populate_fromlocus(locid, gendict[locid])


def populate_anticrispr(data_dir):
    with open(os.path.join(data_dir, 'antiCRISPR_accessions.txt')) as f:
        accession_list = list(read_accession_file(f))
    print("Fetching AntiCRISPR entries")
    result_handle = Entrez.efetch(
        db='protein', rettype="fasta", id=accession_list)
    for record in tqdm(SeqIO.parse(result_handle, 'fasta')):
        spacer, _ = AntiCRISPR.objects.get_or_create(
            accession=record.name,
            sequence=str(record.seq))
        spacer.save()


def main(email, data_dir, limit):
    Entrez.email = email
    print("Starting organism population")
    populate_organism(data_dir)
    print("Starting LSR population")
    populate_lsrpair(data_dir, limit)
    print("Starting AntiCRISPR population")
    populate_anticrispr(data_dir)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Populate the phageParser database with data from NCBI')
    parser.add_argument('email', nargs=1,
                        help='your email address (does not need to be registered, just used to identify you)')
    parser.add_argument('data_dir', nargs=1, default='data',
                        help='a directory for input and temporary files')
    parser.add_argument('test_limit', nargs=1, default=-1,
                        help='useful for limiting the size of the test set')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), args.data_dir))
    main(args.email, data_dir, args.limit)
