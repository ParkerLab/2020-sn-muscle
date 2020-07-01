#!/usr/bin/env python
# coding: utf-8

import pybedtools as bt
import csv
import sys
import argparse
import logging
import pandas as pd

parser = argparse.ArgumentParser(description='Given a BED file with at least 6 columns (representing gene starts and ends), spit out a feature file for liger.', add_help = True)
parser.add_argument('bed', type = str,  help = 'BED file.')
parser.add_argument('chrom_sizes', type = str,  help = 'Chromosome size file.')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')


def list_to_bedtool(l):
    return bt.BedTool('\n'.join(l), from_string = True)


# remove genes that dont appear on chromosomes in the chromosome size file...
chrom_sizes = pd.read_csv(args.chrom_sizes, delimiter='\t', header=None, names=['chrom', 'size'])
chromosomes = set(chrom_sizes.chrom)

logging.info('Checking that each gene is only assigned to one chromosome and strand...')
bed = pd.read_csv(args.bed, delimiter='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
bed_unique_strands = bed[['chrom', 'name', 'strand']].drop_duplicates()
counts = bed_unique_strands['name'].value_counts()
duplicates = counts[counts > 1].index.tolist()
for d in duplicates:
    logging.info('Dropping {} (assigned either to multiple chromosomes or multiple strands)'.format(d))
bed = bed[~bed['name'].isin(duplicates)]

# drop chromosomes that aren't in chrom size file
missing_chroms = set([i for i in bed.chrom if i not in chromosomes])
for m in missing_chroms:
    logging.info('Dropping chromosome {} (not in chromosome size file)'.format(m))
bed = bed[bed.chrom.isin(chromosomes)]

features = {name: [] for name in bed['name']}
strands = {row['name']: row['strand'] for index, row in bed.iterrows()}

logging.info('Splitting by features...')
for index, row in bed.iterrows():
    features[row['name']].append('{chrom} {start} {end} {name} . {strand}'.format(**row))


logging.info('Merging within each feature')

genes = []
count = 0
total_features = len(features)

for feature, regions in features.items():
    count += 1
    if count % 1000 == 0:
        logging.info('Processed {} of {} features'.format(count, total_features))
    f = list_to_bedtool(regions).sort().merge(s=True, c=[4,5,6], o=['distinct','distinct','distinct']).slop(g=args.chrom_sizes, l=3000, r=0, s=True).to_dataframe()
    assert(len(f.columns) == 6)
    f.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    genes.append(f)

genes = pd.concat(genes)

# remove non-contiguous features...
counts = genes['name'].value_counts()
duplicates = counts[counts > 1].index.tolist()
for d in duplicates:
    logging.info('Dropping {} (non-contiguous)'.format(d))
genes = genes[~genes.name.isin(duplicates)]
genes[['chrom', 'start', 'end', 'name']].to_csv(sys.stdout, sep='\t', index=False)

logging.info('Done')
