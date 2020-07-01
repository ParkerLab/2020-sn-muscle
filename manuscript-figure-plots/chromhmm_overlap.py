#!/usr/bin/env python3
import os
import pybedtools as bt
import logging
import argparse

parser = argparse.ArgumentParser(description='Prepare a directory containing everything needed for a GB session', add_help = True)
parser.add_argument('peaks', type = str,  help = 'Bed file of peaks.')
parser.add_argument('tss', type = str,  help = 'Bed file of TSS (not gzipped)')
parser.add_argument('chrom_sizes', type = str,  help = 'Chrom sizes file')
parser.add_argument('chromhmm', nargs = '+', type = str, help = 'Bed files (each only 3 columns) of chromHMM regions. File names should follow format: {tissue}.{state}.bed')
parser.add_argument('--tss-extension', dest = 'tss_extension', type = int, default = 1000, help = 'Distance from TSS (in bp) to be considered TSS-distal/proximal (default: 1000).')

args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s: %(message)s')


def parse_chromhmm_file_name (f):
    """ChromHMM files are assumed to following the scheme: {tissue}.{state}.bed"""
    tissue, state, suffix = os.path.basename(f).split('.')
    return {'tissue': tissue, 'state': state}


logging.info('Extending TSS by {}'.format(args.tss_extension))
tss_extended = bt.BedTool(args.tss).slop(b = args.tss_extension, g = args.chrom_sizes).sort()

logging.info('Creating TSS-proximal and TSS-distal peaks')
proximal = bt.BedTool(args.peaks).sort().intersect(tss_extended, wa = True).merge()
distal = bt.BedTool(args.peaks).sort().intersect(tss_extended, wa = True, v = True).merge()

logging.info('Overlapping distal intervals with chromHMM files')
count = 0
distal_sums = dict()
for f in args.chromhmm:
    count += 1
    logging.info('Overlapping with chromHMM file {} of {} ({})'.format(count, len(args.chromhmm), f))
    distal_overlap = distal.intersect(b = f, wo = True)
    logging.info('Calculating overlap sums')
    for i in distal_overlap:
        chromhmm_file = f
        amount_of_overlap = int(i[6])
        if chromhmm_file not in distal_sums:
            distal_sums[chromhmm_file] = 0
        distal_sums[chromhmm_file] += amount_of_overlap

logging.info('Overlapping proximal intervals with chromHMM files')
count = 0
proximal_sums = dict()
for f in args.chromhmm:
    count += 1
    logging.info('Overlapping with chromHMM file {} of {} ({})'.format(count, len(args.chromhmm), f))
    proximal_overlap = proximal.intersect(b = f, wo = True)
    logging.info('Calculating overlap sums')
    for i in proximal_overlap:
        chromhmm_file = f
        amount_of_overlap = int(i[6])
        if chromhmm_file not in proximal_sums:
            proximal_sums[chromhmm_file] = 0
        proximal_sums[chromhmm_file] += amount_of_overlap

# output
print('tissue\tstate\toverlap\ttss_relative')
for chromhmm_file, overlap in distal_sums.items():
    print('{}\t{}\t{}\tdistal'.format(parse_chromhmm_file_name(chromhmm_file)['tissue'], parse_chromhmm_file_name(chromhmm_file)['state'], overlap))
for chromhmm_file, overlap in proximal_sums.items():
    print('{}\t{}\t{}\tproximal'.format(parse_chromhmm_file_name(chromhmm_file)['tissue'], parse_chromhmm_file_name(chromhmm_file)['state'], overlap))

logging.info('Done')
