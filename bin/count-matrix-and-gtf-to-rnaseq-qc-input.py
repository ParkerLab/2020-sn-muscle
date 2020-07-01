#!/usr/bin/env python

import os
import sys
import re
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

COUNT_MATRIX, GTF = sys.argv[1:]

# map gene name --> chromosome
logging.info('Reading gene -> chrom mappings')
gene_to_chromosome = dict()

with open(GTF, 'r') as f:
    for line in f:
        chrom = line.rstrip().split('\t')[0]
        gene_name = re.search(r'gene_name "(.*?)";', line).group(1)
        gene_to_chromosome[gene_name] = chrom

for gene, chrom in gene_to_chromosome.items():
    logging.info('Mapped gene {} to chrom {}'.format(gene, chrom))

logging.info('Gathering stats')
barcode_to_stats = dict() # barcode --> number_genes, number_umis, number_mitochondrial_umis

count = 0
with open(COUNT_MATRIX, 'r') as f:
    for line in f:
        count += 1
        if count % 1000000 == 0:
            logging.info('Processed {} lines'.format(count))
        library, barcode, gene, count = line.rstrip().split('\t')
        count = int(count)
        if barcode not in barcode_to_stats:
            barcode_to_stats[barcode] = [0, 0, 0]
        barcode_to_stats[barcode][0] += 1
        barcode_to_stats[barcode][1] += count
        if gene_to_chromosome[gene] == 'chrM':
            barcode_to_stats[barcode][2] += count

print('barcode\tnumber_genes\tnumber_umis\tfraction_mitochondrial')

for barcode, stats in barcode_to_stats.items():
    print('{}\t{}\t{}\t{}'.format(barcode, stats[0], stats[1], stats[2] / stats[1]))
