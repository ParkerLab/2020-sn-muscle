#!/usr/bin/env python

import os
import sys
import pybedtools as bt
import pandas as pd

BEDS = sys.argv[1:]

beds = pd.concat([pd.read_csv(f, usecols=[0, 1, 2], delimiter='\t', header=None, names=['chrom', 'start', 'end']) for f in BEDS])
merged = bt.BedTool().from_dataframe(beds).sort().merge()
master_peak_df = merged.to_dataframe().assign(overlaps=0, peak=lambda df: df.chrom + ':' + df.start.astype(str) + ':' + df.end.astype(str)).set_index('peak')
sys.stderr.write('Found {} master peaks\n'.format(len(master_peak_df)))
sys.stderr.flush()

for f in BEDS:
    sys.stderr.write('Finding overlaps with {}\n'.format(f))
    sys.stderr.flush()
    overlapped = merged.intersect(f, u=True).to_dataframe().applymap(str)
    overlapped = overlapped.chrom + ':' + overlapped.start + ':' + overlapped.end
    master_peak_df.loc[overlapped,'overlaps'] = master_peak_df.loc[overlapped,'overlaps'] + 1

master_peak_df[['chrom', 'start', 'end', 'overlaps']].to_csv(sys.stdout, sep='\t', index=False, header=False)
