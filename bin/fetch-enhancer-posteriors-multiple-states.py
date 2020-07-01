#!/usr/bin/env python
import pybedtools as bt
import pandas as pd
import argparse
import os
import logging
import sys
import gzip
import itertools

logging.basicConfig(level=logging.DEBUG, format = '%(asctime)s - %(levelname)s: %(message)s')

parser = argparse.ArgumentParser()
parser.add_argument('bed', help = 'Bed file of regions for which to fetch posteriors. Should contain at least 4 columns.')
parser.add_argument('roadmap_posterior_file', help = 'Roadmap posterior file (e.g., E016_15_coreMarks_chr18_posterior.txt.gz).')
parser.add_argument('state', nargs = '+', help = 'States for with posteriors should be extracted.')
args = parser.parse_args()

# take the max in each master peak as the posterior...
def read_posterior_file_header(posterior_file):
    cell_type = None
    chrom = None
    with gzip.open(posterior_file, 'rt') as f:
        cell_type, chrom = f.readline().rstrip().split()
    return {'chrom': chrom, 'cell_type': cell_type}


def load_posterior_file (posterior_file, states):
    # return a bedtool -- chrom, start, end, ., score
    df = pd.read_csv(posterior_file, delimiter='\t', usecols=states, skiprows=1)
    df['score'] = df.max(axis=1)
    starts = itertools.count(0, 200)
    df['chrom'] = read_posterior_file_header(posterior_file)['chrom']
    df['start'] = [next(starts) for i in range(len(df))]
    df['end'] = df.start + 200
    df['name'] = 'x'
    return bt.BedTool().from_dataframe(df[['chrom', 'start', 'end', 'name', 'score']])

logging.info('Reading bed file...')
bed = bt.BedTool().from_dataframe(pd.read_csv(args.bed, delimiter='\t', header=None, names=['chrom', 'start', 'end', 'name'], usecols=[0, 1, 2, 3]))
logging.info('Finished reading bed file.')

out = []

for posterior_file in [args.roadmap_posterior_file]:
    state_data = load_posterior_file(posterior_file, args.state)
    cell_type = read_posterior_file_header(posterior_file)['cell_type']
    logging.info('Fetching posteriors in bed regions...')
    overlaps = bed.intersect(state_data, wo = True).cut([0, 1, 2, 3, 8]).to_dataframe()
    overlaps['score'] = overlaps['score'].map(float)
    maxes = overlaps.groupby(by=['name']).score.agg(max).reset_index() # name is a peak name
    maxes['cell_type'] = cell_type
    out.append(maxes[['name', 'cell_type', 'score']])
    logging.info('Finished with {}'.format(posterior_file))

logging.info('Printing (unsorted) output...')
pd.concat(out).sort_values(['name', 'cell_type']).to_csv(sys.stdout, sep='\t', index=False)

