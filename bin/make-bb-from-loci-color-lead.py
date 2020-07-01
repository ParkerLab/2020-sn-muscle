#!/usr/bin/env python
import os
import sys
import pandas as pd

locus_file, index_snp_list = sys.argv[1:]
# index SNP list should be format chrom\tpos

index_snps = []
with open(index_snp_list, 'r') as f:
    for line in f:
        chrom, pos = line.rstrip().split()
        pos = int(pos)
        index_snps.append('{}:{}:{}'.format(chrom, pos-1, pos))

df = pd.read_csv(locus_file, delimiter='\t', header=None, names=['chrom', 'start', 'end', 'locus'])
df = df.sort_values(['chrom', 'start'])
df['name'] = df.chrom + ':' + df.start.map(str) + ':' + df.end.map(str)
df['score'] = 1
df['strand'] = '.'
df['count'] = 1
df['sizes'] = df.end - df.start
df['starts'] = 0
df['color'] = df.name.map(lambda x: '128,0,128' if x in index_snps else '255,0,0')
df[['chrom', 'start', 'end', 'name', 'score', 'strand', 'start', 'end', 'color', 'count', 'sizes', 'starts']].to_csv(sys.stdout, sep='\t', header=None, index=False)
