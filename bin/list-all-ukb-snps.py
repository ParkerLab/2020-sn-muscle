#!/usr/bin/env python

import os
import sys
import pandas as pd
import gzip
import glob
import csv
import re
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

conversions = set()
GWAS_FILES = sys.argv[1:]

for i in GWAS_FILES:
    logging.info('Reading file {}'.format(i))
    with gzip.open(i, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for line in reader:
            conversions.add(line['variant'])

for key in conversions:
    print(key)
