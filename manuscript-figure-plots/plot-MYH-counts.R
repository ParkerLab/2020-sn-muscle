#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)

args <- commandArgs(T)
CLUSTERS <- args[1]
COUNTS <- args[2]
CLUSTER_NAMES <- args[3]
LIBRARY_LABELS <- args[4]

library_to_modality_and_species <- read.table(LIBRARY_LABELS, head=T, as.is=T, sep='\t') %>%
  dplyr::select(library, species, modality)

clusters <- read.table(CLUSTERS, head=F, sep='\t', col.names = c('library', 'barcode', 'cluster'), colClasses = c('character'))
counts <- read.table(COUNTS, head=F, sep='\t', col.names = c('library', 'barcode', 'gene', 'count'), colClasses = c('character', 'character', 'character', 'numeric'))
cluster_names <- read.table(CLUSTER_NAMES, head=T, sep='\t', stringsAsFactors = F, colClasses = c('character'))
colnames(cluster_names) <- c('cluster', 'cluster_name')
cluster_names$cluster_name <- factor(cluster_names$cluster_name, levels=rev(cluster_names$cluster_name), ordered=T)
counts <- left_join(clusters, counts) %>%
  left_join(cluster_names)
counts <- left_join(counts, library_to_modality_and_species)
counts <- counts %>%
  tidyr::spread(key=gene, value=count, fill=0)

counts <- counts[counts$modality=='ATAC',]

ALPHAS <- c('Rat'=0.3, 'Human'=0.05)

p <- ggplot(counts[grep('fiber', counts$cluster_name),]) +
  geom_point(aes(x=MYH1+MYH2+MYH4, y=MYH7, alpha=species), stroke=0) +
  theme_bw() +
  facet_grid(species~cluster_name) +
  scale_alpha_manual(values=ALPHAS) +
  guides(alpha=F)
png('MYH1-plus-MYH4-vs-MYH7.png', height=4, width=6, units='in', res=300)
p
dev.off()

p <- ggplot(counts) +
  geom_point(aes(x=MYH1+MYH2+MYH4, y=MYH7), alpha=0.2, stroke=0) +
  theme_bw() +
  facet_grid(cluster_name~species)
png('MYH1-plus-MYH4-vs-MYH7-all-cell-types.png', height=11, width=5, units='in', res=300)
p
dev.off()

