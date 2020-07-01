#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(glue)

args <- commandArgs(T)
DIM_FILE <- args[1]
CLUSTER_FILE <- args[2]
OUT <- args[3]

tmp <- read.table(DIM_FILE, head = F, sep = '\t', col.names = c('library', 'barcode', 'dim1', 'dim2'), colClasses = c('character', 'character', 'numeric', 'numeric'))
tmp <- tmp[sample(1:nrow(tmp)),]

clusters <- read.table(CLUSTER_FILE, head = F, sep = '\t', col.names = c('library', 'barcode', 'cluster'), colClasses = c('character', 'character', 'numeric'))

tmp <- left_join(tmp, clusters) %>%
	dplyr::group_by(cluster) %>%
	dplyr::mutate(size=n()) %>%
	dplyr::ungroup() %>%
	dplyr::select(library, barcode, dim1, dim2, cluster, size)
tmp <- tmp[rev(order(tmp$cluster)),]
tmp$label <- paste0(tmp$cluster, ' (n = ', tmp$size, ')')
tmp$label <- factor(tmp$label, levels=unique(tmp$label), ordered=T)

p <- ggplot(tmp) +
  geom_point(aes(x = dim1, y = dim2, color = as.factor(label)), size = 1) +
  ylab('Dim. 2') +
  xlab('Dim. 1') +
  theme_bw() +
  scale_color_viridis_d() +
  guides(color = guide_legend(title='Cluster'))
pdf(OUT, height = 5, width = 6)
p
dev.off()
