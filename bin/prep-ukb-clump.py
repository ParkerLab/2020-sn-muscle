#!/usr/bin/env python

import os
import sys
import csv

gwas = sys.argv[1]

print('SNP\tP')

with open(gwas, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for line in reader:
        if 'rs' in line['rsID'] and line['pval'] != '':
            print('{}\t{}'.format(line['rsID'], line['pval']))
