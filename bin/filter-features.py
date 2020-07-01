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
parser.add_argument('--keep-features', dest = 'keep_features', default = None, type = str, action='store', help = 'Keep features listed in this file.')
parser.add_argument('--keep-nuclei', dest = 'keep_nuclei', default = None, type = str, action='store', help = 'Keep nuclei listed in this file.')
args = parser.parse_args()

nuclei = set()
if args.keep_nuclei:
    with open(args.keep_nuclei, 'r') as f:
        for line in f:
            library, barcode = line.rstrip().split()
            nuclei.add('{}-{}'.format(library, barcode))

features = set()
if args.keep_features:
    with open(args.keep_features, 'r') as f:
        for line in f:
            feature = line.rstrip()
            features.add(feature)


with open(args.feature_file, 'r') as f:
    for line in f:
        library, barcode, feature, score = line.rstrip().split('\t')
        nucleus = '{library}-{barcode}'.format(**locals())
        if args.keep_features and not feature in features:
            continue
        if args.keep_nuclei and not nucleus in nuclei:
            continue
        print(line.rstrip())
