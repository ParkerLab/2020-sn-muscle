#!/usr/bin/env python

import os
import sys

overlaps, library_name = sys.argv[1:]

counts = dict() # barcode --> feature --> count

with open(overlaps, 'r') as f:
    for line in f:
        barcode, feature = line.rstrip().split('\t')
        if barcode not in counts:
            counts[barcode] = dict()
        if feature not in counts[barcode]:
            counts[barcode][feature] = 0
        counts[barcode][feature] += 1

for barcode in sorted(counts.keys()):
    for feature, count in sorted(counts[barcode].items()):
        print('{}\t{}\t{}\t{}'.format(library_name, barcode, feature, count))
