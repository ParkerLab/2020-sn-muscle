#!/usr/bin/env python

import os
import sys
import json
import pandas as pd

ROOT = sys.argv[1]
SAMPLE_INFO = os.path.join(ROOT, 'sample_info', 'sample_info.txt')

def get_samples(library):
    if 'rn6' in library:
        return ['rat1']
    elif '63_20' in library or '63_40' in library:
        return ['KSM1', 'KSM2']
    else:
        return ['KSM1']


sample_info = pd.read_csv(SAMPLE_INFO, delimiter='\t')
atac = sample_info[sample_info.modality=='ATAC']
rna = sample_info[sample_info.modality=='RNA']
atac_libraries = []
rna_libraries = []

for index, row in atac.iterrows():
    genomes = row['genomes'].split(',')
    for g in genomes:
        atac_libraries.append('{}-{}'.format(row['library'], g))

for index, row in rna.iterrows():
    genomes = row['genomes'].split(',')
    for g in genomes:
        rna_libraries.append('{}-{}'.format(row['library'], g))

CONFIG = {'libraries': dict()}

for lib in atac_libraries:
    genome = lib.split('-')[1]
    CONFIG['libraries'][lib] = {
            'genome': genome, 
            'samples': get_samples(lib), 
            'modality': 'ATAC', 
            # 'counts': os.path.join(ROOT, 'work', 'atacseq', 'results', 'gene-counts', '{}.counts.txt'.format(lib)), 
            'counts': os.path.join(ROOT, 'work', 'counts-atac', 'results', 'uncorrected-counts', '{}.counts.txt'.format(lib)), 
            'ataqv_json': os.path.join(ROOT, 'work', 'atacseq', 'results', 'ataqv', '{}.ataqv.json.gz'.format(lib)), 
            'pruned': os.path.join(ROOT, 'work', 'atacseq', 'results', 'prune', '{}.pruned.bam'.format(lib))}

for lib in rna_libraries:
    genome = lib.split('-')[1]
    CONFIG['libraries'][lib] = {
            'genome': genome, 
            'samples': get_samples(lib), 
            'modality': 'RNA', 
            'starsolo': [os.path.join(ROOT, 'work', 'rnaseq', 'results', 'starsolo', lib, 'Solo.out', 'GeneFull', 'raw', x) for x in ['features.tsv', 'barcodes.tsv', 'matrix.mtx']],
            #'counts': os.path.join(ROOT, 'work', 'compare-cellbender', 'results', 'corrected-counts', lib + '.features.txt'),
            'counts': os.path.join(ROOT, 'work', 'counts', 'results', 'uncorrected-counts', lib + '.features.txt'),
            'qc': os.path.join(ROOT, 'work', 'rnaseq', 'results', 'qc', '{}.qc.json.gz'.format(lib)),
            'pruned': os.path.join(ROOT, 'work/rnaseq/results/prune/{}.before-dedup.bam'.format(lib))}

print(json.dumps(CONFIG, indent=4, sort_keys = True))
