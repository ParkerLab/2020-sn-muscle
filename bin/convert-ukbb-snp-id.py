#!/usr/bin/env python

import os
import sys
import pandas as pd
import gzip

UKBB_GWAS_FILE, CONVERSIONS_FILE, OUT = sys.argv[1:]

conversions = pd.read_csv(CONVERSIONS_FILE, delimiter='\t').set_index('snp')
gwas = None
with gzip.open(UKBB_GWAS_FILE, 'rt') as f:
    gwas = pd.read_csv(f, delimiter='\t')

gwas.n_complete_samples = gwas.n_complete_samples.map(int)
gwas = gwas.rename(columns={'variant': 'snp', 'n_complete_samples': 'N', 'pstat': 'p'}).set_index('snp').join(conversions).reset_index()
snp_info = [x.split(':') for x in gwas.snp]
gwas['A1'] = [x[2] for x in snp_info]
gwas['A2'] = [x[3] for x in snp_info]
gwas['chrom'] = [x[0] for x in snp_info]

gwas = gwas.drop(columns=['snp'])
noconversion = gwas[gwas.rsID.isnull()]
gwas = gwas[~gwas.rsID.isnull()]


gwas.to_csv(OUT, sep='\t', index=False)
noconversion.to_csv(OUT + '.dropped', sep='\t', index=False)
