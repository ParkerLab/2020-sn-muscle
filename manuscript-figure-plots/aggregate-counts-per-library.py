#!/usr/bin/env python

import sys

count_files = sys.argv[1:]

counts = dict() # library -> feature -> count

for c in count_files:
    with open(c, 'r') as f:
        for line in f:
            library, barcode, feature, count = line.rstrip().split()
            count = int(count)
            if library not in counts:
                counts[library] = dict()
            if feature not in counts[library]:
                counts[library][feature] = 0
            counts[library][feature] += count

for library in counts:
    for feature, count in counts[library].items():
        print('{}\t{}\t{}'.format(library, feature, count))
