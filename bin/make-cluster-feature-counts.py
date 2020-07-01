#!/usr/bin/env python

import os
import sys
import logging
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Given a cluster file, sample info file, and a list of feature count files, get per-cluster, genome, modality feature counts and tpms (cluster-genome, modality, feature, count).', add_help = True)
parser.add_argument('clusters', type = str,  help = 'library, barcode, cluster.')
parser.add_argument('sample_info', type = str,  help = 'sample info file, including a library, genomes, and modality columns.')
parser.add_argument('feature_counts', type = str, nargs='+', help = 'List of feature count files')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

library_to_modality = dict()
nucleus_to_cluster = dict()

logging.info('Getting library-to-modality mappings')
for index, row in pd.read_csv(args.sample_info, delimiter='\t').iterrows():
    for genome in row['genomes'].split(','):
        library_to_modality['{}-{}'.format(str(row['library']), genome)] = row['modality']

logging.info('Getting nucleus-to-cluster mappings')
with open(args.clusters, 'r') as f:
    for line in f:
        library, barcode, cluster = line.rstrip().split()
        genome = library.split('-').pop()
        nucleus_to_cluster['{}-{}'.format(library, barcode)] = '{}-{}'.format(cluster, genome)

logging.info('Reading counts')
counts = {cluster: {modality: dict() for modality in library_to_modality.values()} for cluster in nucleus_to_cluster.values()} # cluster --> feature --> count

for feature_counts_file in args.feature_counts:
    logging.info('Reading counts from file {}'.format(feature_counts_file))
    with open(feature_counts_file, 'r') as f:
        for line in f:
            library, barcode, feature, count = line.rstrip().split()
            nucleus = '{}-{}'.format(library, barcode)
            modality = library_to_modality[library]
            count = int(count)
            if nucleus not in nucleus_to_cluster:
                continue
            cluster = nucleus_to_cluster[nucleus]
            counts[cluster][modality][feature] = count if feature not in counts[cluster][modality] else counts[cluster][modality][feature] + count


for cluster, modality_dict in counts.items():
    for modality, feature_dict in modality_dict.items():
        for feature, count in feature_dict.items():
            print('{}\t{}\t{}\t{}'.format(cluster, modality, feature, count))
