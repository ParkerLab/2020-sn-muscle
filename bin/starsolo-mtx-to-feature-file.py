#!/usr/bin/env python

import os
import sys

item1, item2, matrix = sys.argv[1:] # e.g. barcodes.tsv, features.tsv, matrix.mtx

item_1 = dict()
item_2 = dict()

with open(item1, 'r') as f:
    line_count = 0
    for line in f:
        line_count += 1
        line = line.rstrip()
        item_1[line_count] = line

with open(item2, 'r') as f:
    line_count = 0
    for line in f:
        line_count += 1
        line = line.rstrip()
        item_2[line_count] = line


with open(matrix, 'r') as f:
    line_count = 0 # current line, excluding header
    for line in f:
        line = line.rstrip()
        if line_count == 0 and '%' in line:
            # header line...
            continue
        line_count += 1
        if line_count == 1:
            # want to skip the first line after the header
            continue
        index_1, index_2, value = [int(x) for x in line.split()]
        print('{}\t{}\t{}'.format(item_2[index_2], item_1[index_1], value))
