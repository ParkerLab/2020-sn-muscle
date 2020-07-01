#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(glue)

args <- commandArgs(T)
DIM_FILE <- args[1]
PREFIX <- args[2]

tmp <- read.table(DIM_FILE, head = F, sep = '\t', col.names = c('library', 'barcode', 'dim1', 'dim2'), colClasses = c('character', 'character', 'numeric', 'numeric'))
tmp <- tmp[sample(1:nrow(tmp)),]

RNA <- paste(c('133155', '133156', '133157', '133158', '63_20_rna', '63_40_rna'), 'hg19', sep='-')

tmp$modality <- ifelse(tmp$library %in% RNA, 'RNA', 'ATAC')

number_libraries <- length(unique(tmp$library))

p <- ggplot(tmp) +
  geom_point(aes(x = dim1, y = dim2, color=modality), alpha=0.1) +
  ylab('Dim. 2') +
  xlab('Dim. 1') +
  theme_bw() +
  scale_color_viridis_d()
pdf(glue('{PREFIX}umap-by-modality.pdf'), height = 5, width=6)
p
dev.off()

tmp$library <- paste('library', tmp$library)

p <- ggplot(tmp) +
  geom_point(aes(x = dim1, y = dim2, color=modality), alpha=0.1) +
  ylab('Dim. 2') +
  xlab('Dim. 1') +
  theme_bw() +
  scale_color_viridis_d() +
  facet_wrap(~library)
pdf(glue('{PREFIX}umap-by-library.pdf'), height = 2*sqrt(number_libraries), width = 1+2*sqrt(number_libraries))
p
dev.off()
