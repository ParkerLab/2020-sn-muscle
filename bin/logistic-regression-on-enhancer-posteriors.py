#!/usr/bin/env python
import argparse
import os
import logging
import sys
import gzip
import pandas as pd
import statsmodels.discrete.discrete_model as sm
from statsmodels.api import add_constant
import numpy as np
import copy
import pybedtools # so that we can subset to distal peaks only...

logging.basicConfig(level=logging.DEBUG, format = '%(asctime)s - %(levelname)s: %(message)s')

parser = argparse.ArgumentParser()
parser.add_argument('accessibility', help = 'Bed file of master peaks (peak name, cluster that the peak is accessible in). If a peak is accessible in multiple clusters, then there will be multiple lines for that peak.')
parser.add_argument('posterior_files', nargs = '+', help = 'List of text files, one per roadmap cell type, giving the posterior for being an enhancer in that master peak in that cell type. (peak, cell type, score)')
args = parser.parse_args()

#accessibility_file = '/lab/work/porchard/PL2689/work/process-by-cluster/results/master-peaks/accessibility.bed'
#posterior_dir = '/lab/work/porchard/PL2689/work/process-by-cluster/results/master-peak-posteriors'
#posterior_files = [os.path.join(posterior_dir, x) for x in os.listdir(posterior_dir) if '.posteriors.bed' in x]

# Determine which master peaks are accessible in each cluster
accessible = dict()
clusters = set()

logging.info('Reading master peak accessibility file')
with open(args.accessibility, 'r') as f:
    for line in f:
        peak, cluster = line.rstrip().split('\t')
        if peak not in accessible:
            accessible[peak] = set()
        accessible[peak].add(cluster)
        clusters.add(cluster)

posteriors = dict() # peak --> cell_type --> posterior
cell_types = set()

for posterior_file in args.posterior_files:
    logging.info('Reading posterior file {}'.format(posterior_file))
    with open(posterior_file, 'r') as f:
        for line in f:
            peak, cell_type, posterior = line.rstrip().split()
            if peak not in accessible:
                continue
            if peak not in posteriors:
                posteriors[peak] = dict()
            posteriors[peak][cell_type] = float(posterior)
            cell_types.add(cell_type)

# Some peaks in the accessibility file may not have scores (e.g., because it's in a non-human species and the peak couldn't be lifted over)
# remove such cases...
NO_SCORE = set()
for peak in accessible:
    if peak not in posteriors:
        NO_SCORE.add(peak)
logging.info('{} of {} peaks are missing scores; dropping these'.format(len(NO_SCORE), len(posteriors)))
for peak in NO_SCORE:
    del accessible[peak]


# for each cluster, model the the accessibility of each peak as a function of the posteriors
models = dict() # cluster --> cell_type --> model

number_models = len(clusters) * len(cell_types)
count = 0
for cluster in clusters:
    models[cluster] = dict()
    for cell_type in cell_types:
        count += 1
        logging.info('Running model {} of {} (cluster {} and cell type {})'.format(count, number_models, cluster, cell_type))
        peaks = accessible.keys()
        accessibility = np.array([1 if cluster in accessible[peak] else 0 for peak in peaks])
        posterior = np.array([posteriors[peak][cell_type] for peak in peaks])
        mod = sm.Logit(accessibility, add_constant(posterior, prepend = False)).fit()
        models[cluster][cell_type] = copy.deepcopy(mod)

# output all the model results
print('cluster\tcell_type\tpvalue\tcoef')
for cluster in clusters:
    for cell_type in cell_types:
        p = models[cluster][cell_type].pvalues[0]
        coef = models[cluster][cell_type].params[0]
        print('{cluster}\t{cell_type}\t{p}\t{coef}'.format(**locals()))
