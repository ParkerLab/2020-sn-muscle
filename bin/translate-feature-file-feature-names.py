#!/usr/bin/env python

import os
import sys
import csv
import argparse
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

parser = argparse.ArgumentParser()
parser.add_argument('--drop-if-missing', dest='drop_if_missing', action='store_true', default = False, help = 'If a feature cannot be translated, drop it.')
parser.add_argument('--keep-if-missing', dest='keep_if_missing', action='store_true', default = False,  help = 'If a feature cannot be translated, keep it.')
parser.add_argument('feature_file', type = str, help = 'The feature file to be translated.')
parser.add_argument('translations', type = str, help = 'The file of feature translations.')
parser.add_argument('--from', dest='fr', type = str, help = 'The name of the translations file corresponding to the old feature names')
parser.add_argument('--to', type = str, help = 'The name of the translations file corresponding to the new feature names')
args = parser.parse_args()

if (not args.drop_if_missing and not args.keep_if_missing) or (args.drop_if_missing and args.keep_if_missing):
    raise ValueError('Must set either the --drop-if-missing or the --keep-if-missing argument')

# read in the translations...
translations = dict()

with open(args.translations, 'r') as f:
    r = csv.DictReader(f, delimiter='\t')
    for line in r:
        fr = line[args.fr]
        to = line[args.to]
        if fr in translations:
            raise ValueError('There are multiple translations from {args.fr} {fr} to {args.to}'.format(**locals()))
        translations[fr] = to

dropped_features = set()
not_translated = set()

with open(args.feature_file, 'r') as f:
    for line in f:
        library, barcode, feature, value = line.rstrip().split('\t')
        
        if feature not in translations and args.drop_if_missing:
            dropped_features.add(feature)
            continue
        elif feature not in translations and args.keep_if_missing:
            not_translated.add(feature)
        else:
            feature = translations[feature]
        
        print('\t'.join([library, barcode, feature, value]))

for feature in dropped_features:
    logging.warning('Dropped feature {feature} (no translation)'.format(**locals()))

for feature in not_translated:
    logging.warning('Keeping feature {feature} as is (no translation)'.format(**locals()))
