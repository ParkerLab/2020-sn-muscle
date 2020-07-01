#!/usr/bin/env python

import logging
import argparse

parser = argparse.ArgumentParser(description='Calculate TSS enrichment from a bam file.', add_help = True)
parser.add_argument('file_to_filter', type = str,  help = 'Text file to filter. Must be tab-separated.')
parser.add_argument('keep', type = str, help = 'Keep lines that have this ID listed in this file.')
parser.add_argument('filter_on_index', type = int,  help = 'Index of column to filter on (indexed from 1).')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

keep = set()
keep_count = 0
drop_count = 0

logging.info('Reading file of items to keep')

with open(args.keep, 'r') as f:
    for line in f:
        line = line.rstrip()
        keep.add(line)

logging.info('Filtering file')

with open(args.file_to_filter, 'r') as f:
    for line in f:
        line = line.rstrip().split('\t')
        if line[args.filter_on_index - 1] in keep:
            print('\t'.join(line))
            keep_count += 1
        else:
            drop_count += 1

logging.info('Finished filtering file')
logging.info('Kept {} and dropped {} items'.format(keep_count, drop_count))
