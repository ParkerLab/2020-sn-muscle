#!/usr/bin/env python
import os
import sys
import pandas as pd

sheet_1 = pd.read_excel('41598_2019_57110_MOESM5_ESM-1.xlsx', skiprows=1, nrows=40, usecols=list(range(4)))
sheet_1.columns = ['gene', 'fiber_type', 'log2FC', 'padj']
sheet_1.to_csv(sys.stdout, sep='\t', index=False)
