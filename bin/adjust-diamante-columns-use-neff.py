#!/usr/bin/env python

import os
import sys
import pandas as pd

SUMSTATS_FILE = sys.argv[1]

gwas = pd.read_csv(SUMSTATS_FILE, delimiter='\t')
gwas = gwas.drop(columns=['SNP']).rename(columns={'Neff': 'N'})
gwas.to_csv(sys.stdout, sep='\t', index=False)
