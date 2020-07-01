#!/usr/bin/env python

import sys
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

items = set()

fs = sys.argv[1:]
for infile in fs:
    logging.info('Processing file {}'.format(infile))
    with open(infile, 'r') as f:
        for line in f:
            line = line.rstrip()
            items.add(line)

for i in sorted(list(items)):
    print(i)

logging.info('Done')

