#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd
import numpy as np

SOUPORCELL_VCF, GENOTYPES, SOUPORCELL_CLUSTERS, LIBRARY, OUT = sys.argv[1:]

# for testing
# SOUPORCELL_VCF = '/lab/work/porchard/sn-muscle-project/work/souporcell/results/sos/63_20_rna/cluster_genotypes.vcf'
# GENOTYPES = '/lab/work/porchard/sn-muscle-project/work/snp-calling-on-bulk/results/genotypes/snp-calls.vcf.gz'
# SOUPORCELL_ASSIGNMENTS = '/lab/work/porchard/sn-muscle-project/work/souporcell/results/sos/63_20_rna/clusters.tsv'
# LIBRARY = '63_20_rna-hg19'

# read the genotypes
genotypes = pd.read_csv(GENOTYPES, delimiter='\t', comment='#', header=None)
genotypes.columns = ['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'KSM1', 'KSM2']
genotypes = genotypes[['chrom', 'pos', 'ref', 'alt', 'KSM1', 'KSM2']]
genotypes.KSM1 = genotypes.KSM1.map(lambda x: x.split(':')[0])
genotypes.KSM2 = genotypes.KSM2.map(lambda x: x.split(':')[0])
genotypes['snp'] = genotypes.chrom + ':' + genotypes.pos.map(str) + ':' + genotypes.ref + ':' + genotypes.alt

souporcell = pd.read_csv(SOUPORCELL_VCF, delimiter='\t', comment='#', header=None)
souporcell.columns = ['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'ind0', 'ind1']
souporcell = souporcell[souporcell['filter'] == '.']
souporcell['ind0_total_counts'] = souporcell.ind0.map(lambda x: sum([int(i) for i in x.split(':')[1:3]]))
souporcell['ind1_total_counts'] = souporcell.ind1.map(lambda x: sum([int(i) for i in x.split(':')[1:3]]))
souporcell.ind0 = souporcell.ind0.map(lambda x: x.split(':')[0])
souporcell.ind1 = souporcell.ind1.map(lambda x: x.split(':')[0])
souporcell = souporcell[souporcell.ind0_total_counts>=10]
souporcell = souporcell[souporcell.ind1_total_counts>=10]
souporcell['snp'] = souporcell.chrom + ':' + souporcell.pos.map(str) + ':' + souporcell.ref + ':' + souporcell.alt

compare = genotypes[['KSM1', 'KSM2', 'snp']].set_index('snp').join(souporcell[['ind0', 'ind1', 'snp']].set_index('snp'))
compare = compare[~compare.ind0.isnull()]
compare = compare[~compare.ind1.isnull()]

concordances = []
for sos_ind in [i for i in compare.columns if 'ind' in i]:
    for ksm in [i for i in compare.columns if 'KSM' in i]:
        matched = [int(i) for i in compare[sos_ind] == compare[ksm]]
        concordances.append([sos_ind, ksm, sum(matched) / len(matched)])
concordances = pd.DataFrame(concordances, columns=['individual', 'KSM', 'concordance']).pivot(index='KSM', columns='individual', values='concordance')

print(concordances)

conversions = {}
for i in concordances.columns:
    conversions[i] = concordances.index[concordances[i]==concordances[i].max()].values[0]
conversions

assignments = pd.read_csv(SOUPORCELL_CLUSTERS, delimiter='\t')
assignments['library'] = LIBRARY
assignments.assignment = 'ind' + assignments.assignment.map(str)
assignments.assignment = assignments.assignment.map(lambda x: conversions[x] if x in conversions else 'doublet')
assignments = assignments[['library', 'barcode', 'assignment']]
assignments.to_csv(OUT, sep='\t', index=False, header=False)
