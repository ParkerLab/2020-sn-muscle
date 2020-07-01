#!/usr/bin/env python

import os
import sys
import pandas as pd

ANNOT_FILE = sys.argv[1] # '/lab/work/porchard/sn-muscle-project/data/ldsc-data/baseline/baselineLD.21.annot.gz'

df = pd.read_csv(ANNOT_FILE, delimiter='\t', compression='gzip')
DROP = [i for i in df.columns if 'Trynka' in i or 'Enhancer' in i or 'Hoffman' in i or 'ENCODE' in i or 'H3K27ac' in i or 'BivFlnk' in i or 'BLUEPRINT' in i or 'Villar' in i or 'Ancient_Sequence_Age_Human_Promoter' in i]
DROP.append('GTEx_eQTL_MaxCPP')
df.drop(columns=DROP).to_csv(sys.stdout, index=False, sep='\t')
