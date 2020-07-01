#!/usr/bin/env python 
# coding: utf-8

import logging
import os
import sys
import re
import argparse

logging.basicConfig(level=logging.DEBUG, format = '%(asctime)s - %(levelname)s: %(message)s')

parser = argparse.ArgumentParser()
parser.add_argument('feature_file', type = str, help = 'Library, barcode, feature, value.')
args = parser.parse_args()

FEATURE_COUNT_MATRIX = args.feature_file


logging.info('Loading counts')
counts = {} # nucleus -> feature -> count

with open(FEATURE_COUNT_MATRIX, 'r') as f:
    for line in f:
        library, barcode, feature, count = line.rstrip().split('\t')
        nucleus = '{library}---{barcode}'.format(**locals())
        count = int(count)
        if nucleus not in counts:
            counts[nucleus] = dict()
        if feature not in counts[nucleus]:
            counts[nucleus][feature] = 0
        counts[nucleus][feature] += count

logging.info('Converting to TPM')

for nucleus, feature_counts in sorted(counts.items()):
    library, barcode = nucleus.split('---')
    total_counts = sum(list(feature_counts.values()))
    for feature, count in sorted(feature_counts.items()):
        tpm = 1e6 * count / total_counts
        print('{library}\t{barcode}\t{feature}\t{tpm}'.format(**locals()))
