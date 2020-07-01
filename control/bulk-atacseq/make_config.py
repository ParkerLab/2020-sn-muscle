#!/usr/bin/env python

import os
import sys
import glob
import re
import pandas as pd
import json

FASTQ_FILES = glob.glob('/lab/data/seqcore/NovaA-1*6/320-NM/Sample*/*.fastq.gz')

def parse_fastq_file_path(f):
    RE = r'^(320-NM-\d+)_.*_(S\d+)_R([12])_001.fastq.gz'
    library, readgroup, read = re.search(RE, os.path.basename(f)).groups()
    return {'library': library, 'readgroup': readgroup, 'read': read, 'file': f}

files = pd.DataFrame([parse_fastq_file_path(f) for f in FASTQ_FILES])
files = files[files.library.isin(['320-NM-1', '320-NM-2', '320-NM-3', '320-NM-4'])]
CONFIG = {'libraries': {library: {'genome': 'hg19', 'readgroups': dict()} for library in files.library.unique()}}

for index, row in files.iterrows():
    if row['readgroup'] not in CONFIG['libraries'][row['library']]['readgroups']:
        CONFIG['libraries'][row['library']]['readgroups'][row['readgroup']] = dict()
    CONFIG['libraries'][row['library']]['readgroups'][row['readgroup']][row['read']] = row['file']

print(json.dumps(CONFIG, sort_keys = True, indent = 4))
