#!/usr/bin/env python
from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
from pybedtools import BedTool
import gzip
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bed-files', nargs='+', type=str, help='the UCSC bed file with the regions that make up your annotation')
    parser.add_argument('--bimfile', type=str, help='plink bim file for the dataset you will use to compute LD scores.')
    parser.add_argument('--annot-file', type=str, help='the name of the annot file to output.')
    parser.add_argument('--tss', type=str, default='', help='Bed file of TSS.')
    parser.add_argument('--chrom-sizes', type=str, default='', help='Chrom size file for genome')

    args = parser.parse_args()

    df_bim = pd.read_csv(args.bimfile,
            delim_whitespace=True, usecols = [0,1,2,3], names = ['CHR','SNP','CM','BP'])
    iter_bim = [['chr'+str(x1), x2, x2] for (x1, x2) in np.array(df_bim[['CHR', 'BP']])]
    bimbed = BedTool(iter_bim)
    df_annot = df_bim
    for bed in args.bed_files:
        nm = os.path.basename(bed).replace('.bed', '')
        annotbed = bimbed.intersect(BedTool(bed).sort())
        bp = [x.start for x in annotbed]
        df_int = pd.DataFrame({'BP': bp, nm: 1})
        df_annot = pd.merge(df_annot, df_int, how='left', on='BP')
    if args.tss != '':
        tss = BedTool(args.tss).cut([0, 1, 2]).sort().slop(b=2000, g=args.chrom_sizes)
        for bed in args.bed_files:
            nm = os.path.basename(bed).replace('.bed', '_tssDistal')
            tss_distal_peaks = BedTool(bed).sort().intersect(tss, v=True)
            annotbed = bimbed.intersect(tss_distal_peaks)
            bp = [x.start for x in annotbed]
            df_int = pd.DataFrame({'BP': bp, nm: 1})
            df_annot = pd.merge(df_annot, df_int, how='left', on='BP')
        
    cols = [os.path.basename(f).replace('.bed', '') for f in args.bed_files]
    if args.tss != '':
        cols = cols + [i + '_tssDistal' for i in cols]
    if args.annot_file.endswith('.gz'):
        with gzip.open(args.annot_file, 'wb') as f:
            df_annot[cols].fillna(0).astype(int).to_csv(f, sep="\t", index=False)
    else:
        df_annot[cols].fillna(0).astype(int).to_csv(args.annot_file, sep="\t", index=False)

