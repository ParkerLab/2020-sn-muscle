#!/usr/bin/env python

import os
import sys
import json
import csv
import re
import glob
import pandas as pd

ROOT, SAMPLE_INFO = sys.argv[1:]

DATA = os.path.join(ROOT, 'data')
FASTQ_DIR = os.path.join(DATA, 'fastq')
FASTQ_FILES = glob.glob(os.path.join(FASTQ_DIR, '*.fastq.gz'))

sample_info = pd.read_csv(SAMPLE_INFO, delimiter='\t')
rna = sample_info[sample_info.modality=='RNA']

def parse_fastq_name(f):
    RE = '^(.*)_([a-z]).(\d+).fastq.gz'
    library, readgroup, read = re.match(RE, os.path.basename(f)).groups()
    return {'library': library, 'readgroup': readgroup, 'read': read}


libraries = set([parse_fastq_name(f)['library'] for f in FASTQ_FILES])
libraries = {library: {'genome': '', 'readgroups': {}} for library in libraries}

for f in FASTQ_FILES:
    info = parse_fastq_name(f)
    library = info['library']
    readgroup = '{library}_{readgroup}'.format(**info)
    if readgroup not in libraries[library]['readgroups']:
        libraries[library]['readgroups'][readgroup] = dict()
    read = info['read']
    libraries[library]['readgroups'][readgroup][read] = f

for row in rna.itertuples():
    libraries[str(row.library)]['genome'] = row.genomes.split(',')

remove = [key for key in libraries if key not in rna.library.map(str).values]
for r in remove:
    del libraries[r]


CONFIG = {
    'libraries': libraries
}
print(json.dumps(CONFIG, sort_keys = True, indent = 4))
