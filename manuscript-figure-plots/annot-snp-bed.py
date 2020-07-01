#!/usr/bin/env python
from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
import pybedtools as bt
import gzip
import os
import sys

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--bed-files', nargs='+', type=str, help='the UCSC bed file with the regions that make up your annotation')
    parser.add_argument('--snp-file', type=str, help='File of SNPs (chrom, start, end).')
    parser.add_argument('--pos-file', type=str, help='File of SNPs (chrom, pos).')
    #parser.add_argument('--has-header', type=bool, default=False, help='SNP bed has header.')
    args = parser.parse_args()

    snp_bed = pd.read_csv(args.snp_file, delimiter='\t', header=None) if args.snp_file else pd.read_csv(args.pos_file, delimiter='\t', header=None)
    keep = [0, 1] if args.pos_file else [0, 1, 2]
    snps = snp_bed[keep]
    if args.pos_file:
        snps.columns = ['chrom', 'end']
        snps.end = snps.end.astype(int)
        snps['start'] = snps.end - 1
    elif args.snp_file:
        snps.columns = ['chrom', 'start', 'end']
        snps.start = snps.start.astype(int)
        snps.end = snps.end.astype(int)
    snps['snp'] = snps.chrom.astype(str) + ':' + snps.end.astype(str)
    snps = snps.set_index('snp')

    for bed in args.bed_files:
        sys.stderr.write('Adding annotations from {}\n'.format(bed))
        nm = os.path.basename(bed).replace('.bed', '')
        overlaps = bt.BedTool().from_dataframe(snps[['chrom', 'start', 'end']]).intersect(bt.BedTool(bed), u=True)
        if overlaps.count() == 0:
            snp_bed[nm] = 0
            continue
        overlaps = overlaps.to_dataframe()
        snps_overlapping = overlaps.chrom.astype(str) + ':' + overlaps.end.astype(str)
        overlap_bool = [int(i) for i in snps.index.isin(snps_overlapping)]
        snp_bed[nm] = overlap_bool
    
    snp_bed.to_csv(sys.stdout, sep="\t", index=False)

