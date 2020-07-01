#!/usr/bin/env python

import sys
import logging
import gzip

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

items = set()

fs = sys.argv[1:]
for infile in fs:
    logging.info('Processing file {}'.format(infile))
    if '.gz' in infile:
        with gzip.open(infile, 'rt') as f:
            for line in f:
                line = line.rstrip()
                items.add(line)
    else:
        with open(infile, 'r') as f:
            for line in f:
                line = line.rstrip()
                items.add(line)

for i in sorted(list(items)):
    print(i)

logging.info('Done')

