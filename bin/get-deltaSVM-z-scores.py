#!/usr/bin/env python
from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
import gzip
import os
import sys
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--deltaSVM-file', type=str, help='The deltaSVM file (chrom:pos:allele effect)')
    parser.add_argument('--snp-file', type=str, help='File of SNPs (chrom, start, end).')
    parser.add_argument('--pos-file', type=str, help='File of SNPs (chrom, pos).')
    args = parser.parse_args()

    # append z-score. ignore direction of effect for now
    logging.info('Reading SNPs of interest')
    snps = set()

    snp_file = args.snp_file if args.snp_file else args.pos_file
    with open(snp_file, 'r') as f:
        for line in f:
            line = line.rstrip().split()
            chrom = line[0]
            pos = line[1] if args.pos_file else line[2]
            snp = '{}:{}'.format(chrom, pos)
            snps.add(snp)

    logging.info('Read {} SNPs'.format(len(snps)))
    all_svm_scores = []
    snp_scores = dict()
    logging.info('Reading deltaSVM scores')

    with open(args.deltaSVM_file, 'r') as f:
        count = 0
        for line in f:
            count += 1
            if count % 1000000 == 0:
                logging.info('Read {} scores so far'.format(count))
            snp, score = line.rstrip().split()
            chrom, pos, allele = snp.split(':')
            chrom = 'chr' + chrom if 'chr' not in chrom else chrom
            snp = '{}:{}'.format(chrom, pos)
            score = float(score)
            all_svm_scores.append(score)
            if snp in snps:
                if snp not in snp_scores:
                    snp_scores[snp] = list()
                snp_scores[snp].append(score)
    all_svm_scores = np.array(all_svm_scores)
    all_svm_score_mean = np.mean(all_svm_scores)
    all_svm_score_sd = np.std(all_svm_scores)
    logging.info('Mean and stddev of all SVM scores are {} and {}'.format(all_svm_score_mean, all_svm_score_sd))

    logging.info('Finished reading deltaSVM scores')
    present_count = 0
    missing_count = 0
    with open(snp_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            line_split = line.split()
            chrom = line_split[0]
            pos = line_split[1] if args.pos_file else line_split[2]
            snp = '{}:{}'.format(chrom, pos)
            if snp in snp_scores:
                present_count += 1
                if len(snp_scores[snp]) > 1:
                    logging.warning('Encountered SNP of interest {} multiple times; perhaps there are multiple alleles there. Keeping most extreme score.'.format(snp))
                z_scores = [(i - all_svm_score_mean) / all_svm_score_sd for i in snp_scores[snp]]
                z_score = None
                for i in z_scores:
                    if z_score is None or abs(i) > abs(z_score):
                        z_score = i
                print(line + '\t' + str(z_score))
            else:
                missing_count += 1
                print(line + '\tNA')
    logging.info('Found deltaSVM scores for {} of {} SNPs'.format(present_count, present_count + missing_count))
    logging.info('Done')
