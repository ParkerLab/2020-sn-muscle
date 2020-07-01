#!/usr/bin/env python

import sys

explained_file = sys.argv[1]
nucleotides = ['A', 'C', 'G', 'T']

with open(explained_file, 'r') as f:
    for line in f:
        snp, gkmexplain_score, scores = line.rstrip().split('\t')
        scores = [i.split(',') for i in scores.split(';')]
        for pos in range(len(scores)):
            for nuc, score in zip(nucleotides, scores[pos]):
                print('{snp}\t{pos}\t{nuc}\t{score}'.format(**locals()))
