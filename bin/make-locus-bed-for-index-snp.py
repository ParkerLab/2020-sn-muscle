#!/usr/bin/env python

import os
import sys
import pandas as pd

ld, index_snp = sys.argv[1:]

df = pd.read_csv(ld, delimiter='\t', header = None)
df.columns = ['chrom_1', 'pos_1', 'snp_1', 'chrom_2', 'pos_2', 'snp_2', 'ld']
df = df[df.ld >= 0.8]
df = df[(df.snp_1 == index_snp) | (df.snp_2 == index_snp)]
df = df[df.snp_1 != df.snp_2]

chromosomes = df.chrom_1.append(df.chrom_2)
pos = df.pos_1.append(df.pos_2)

bed = pd.DataFrame({'chrom': chromosomes, 'start': pos - 1, 'end': pos})
bed = bed.sort_values('start').drop_duplicates()
bed.chrom = bed.chrom.map(lambda x: 'chr' + str(x) if 'chr' not in str(x) else x)
assert(bed.chrom.nunique() == 1)
locus = '{}:{}:{}'.format(bed.chrom.values[0], bed.start.min(), bed.end.max())
bed['locus'] = locus

bed.to_csv(sys.stdout, sep='\t', index=False, header=None)
