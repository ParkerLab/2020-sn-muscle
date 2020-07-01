#!/usr/bin/env python

import os
import re
import sys
import gzip
import pandas as pd
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')


gtf, assembly_report = sys.argv[1:]

# read in the chromosome names...
chromosomes = pd.read_csv(assembly_report, delimiter='\t', comment='#', header=None)
chromosomes.columns = ['sequence_name', 'sequence_role', 'assigned_molecule', 'assigned_molecule_location_type', 'genbank', 'relationship', 'refseq', 'assembly_unit', 'sequence_length', 'ucsc']
chromosome_translations = {row['refseq']: row['ucsc'] for index, row in chromosomes.iterrows()}

# for each item in the GTF:
# convert chromosome name
# check the source. if it's Curated Genomic or BestRefSeq, keep
# keep genes 
# as the name of the feature -- use 'gene' attribute

HEADER = re.compile("^#")

with gzip.open(gtf, 'rt') as f:
    for line in f:
        if HEADER.match(line):
            continue
        chrom, source, feature_type, start, end, score, strand, phase, info = line.rstrip().split('\t')
        if feature_type != 'gene':
            continue
        if source != 'BestRefSeq' and source != 'BestRefSeq%2CGnomon' and source != 'Curated Genomic':
            continue
        chrom = chromosome_translations[chrom]
        info = [i.split(' ') for i in info.rstrip(';').split('"; ')]
        info = {i[0]: i[1].replace('"', '') for i in info}
        if info['gene_biotype'] not in ['protein_coding', 'lncRNA']:
            continue
        gene_name = info['gene']
        if chrom == 'na':
            logging.info('Skipping {} (no chromosome translation)'.format(gene_name))
            continue
        print('{chrom}\t{start}\t{end}\t{gene_name}\t.\t{strand}'.format(**locals()))


