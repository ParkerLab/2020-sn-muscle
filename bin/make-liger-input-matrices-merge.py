#!/usr/bin/env python

import os
import sys
import pandas as pd

LIBRARY_LABELS_FILE, FEATURE_FILE, INDIVIDUALS_FILE = sys.argv[1:]

individuals = pd.read_csv(INDIVIDUALS_FILE, delimiter='\t', header=None)
individuals.columns = ['library', 'barcode', 'individual']
individuals['nucleus'] = individuals.library + '-' + individuals.barcode
individuals = individuals.set_index('nucleus').drop(columns=['library', 'barcode'])

library_labels = pd.read_csv(LIBRARY_LABELS_FILE, delimiter='\t').set_index('library')[['modality']]

counts = pd.read_csv(FEATURE_FILE, delimiter='\t', header=None)
counts.columns = ['library', 'barcode', 'gene', 'count']
counts = counts[~counts.gene.isnull()]
counts['nucleus'] = counts.library + '-' + counts.barcode
counts = counts.set_index('nucleus').join(individuals).reset_index().set_index('library').join(library_labels).reset_index()

groups = counts.groupby(['modality', 'individual'])

for name, group in groups:
    fname = '-'.join(name).replace(' ', '') + '.liger.txt'
    group[['nucleus', 'gene', 'count']].pivot(index='gene', columns='nucleus', values='count').fillna(0).applymap(int).to_csv(fname, sep='\t')
