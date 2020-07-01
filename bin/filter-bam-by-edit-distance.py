#!/usr/bin/env python2

import sys
import copy
import os
import json
import gzip
import numpy
import argparse
import re
import pysam
import random
import logging

parser = argparse.ArgumentParser(description='Filter a bam file by read edit distance (must have NM tag)')
parser.add_argument('bam_in', help = 'Bam file to adjust tags in')
parser.add_argument('bam_out', help = 'Output bam')
parser.add_argument('--max-edit-distance', dest = 'max_edit_distance', type = int, default = 99999, help = 'Maximum edit distance to keep.')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

# output new file 
logging.info(f'Writing new bam file ({args.bam_out})')
with pysam.AlignmentFile(args.bam_in, 'rb') as old_bam:
    with pysam.AlignmentFile(args.bam_out, 'wb', template=old_bam) as new_bam:
        keep_count = 0
        remove_count = 0
        for read in old_bam.fetch(until_eof = True):
            if (keep_count + remove_count) % 1000000 == 0:
                logging.info('Processed {} reads...'.format(keep_count + remove_count))
            if read.get_tag('NM') <= args.max_edit_distance:
                keep_count += 1
                new_bam.write(read)
            else:
                remove_count += 1

logging.info(f'Finished filtering {args.bam_in}. Kept {keep_count} reads and removed {remove_count}')
