#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(glue)

args <- commandArgs(T)
UMAP <- args[1]
CLUSTERS <- args[2]

umap <- read.table(UMAP, head = F, sep = '\t', col.names = c('library', 'barcode', 'dim1', 'dim2'), colClasses = c('character', 'character', 'numeric', 'numeric'))
clusters <- read.table(CLUSTERS, head = F, sep = '\t', col.names = c('library', 'barcode', 'cluster'), colClasses = c('character', 'character', 'character'))

both <- left_join(umap, clusters) %>%
  dplyr::group_by(cluster) %>%
  #dplyr::mutate(cluster_size=n(),
  dplyr::mutate(cluster_size=n(),
                cluster_size=formatC(cluster_size, format="d", big.mark=","),
                cluster_label=glue('{cluster} (n={cluster_size})')) %>%
  dplyr::ungroup()
both <- both[order(both$cluster),]
both$cluster_label <- factor(both$cluster_label, levels = unique(both$cluster_label), ordered = T)
both <- both[sample(1:nrow(both)),]

p <- ggplot(both) +
  geom_point(aes(x = dim1, y = dim2, color=cluster_label), alpha=0.3, size=0.1) +
  theme_bw() +
  xlab('UMAP dim. 1') +
  ylab('UMAP dim. 2') +
  guides(color=guide_legend(title='Cluster'))
png('full-umap-legend.png', height=5, width=7, res=300, units='in')
p
dev.off()

p <- ggplot(both) +
  geom_point(aes(x = dim1, y = dim2, color=cluster_label), alpha=0.3, size=0.1) +
  theme_bw() +
  xlab('UMAP dim. 1') +
  ylab('UMAP dim. 2') +
  guides(color=F)
png('full-umap.png', height=3, width=3, res=300, units='in')
p
dev.off()

both$species <- ifelse(gsub('.*-', '', both$library)=='hg19', 'human', 'rat')
both$modality <- ifelse(gsub('-.*', '', both$library) %in% c('125589', '133151', '133152', '133153', '133154', '63_20', '63_40'), 'ATAC', 'RNA') 
both$type <- paste(both$species, both$modality)
both <- both %>%
  dplyr::group_by(type) %>%
  dplyr::mutate(count=n(),
                count=formatC(count, format="d", big.mark=',')) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(type=glue('{type}\n(n={count})'))

p <- ggplot(both) +
  geom_point(aes(x = dim1, y = dim2, color=cluster_label), alpha=0.3, size=0.2) +
  theme_bw() +
  xlab('UMAP dim. 1') +
  ylab('UMAP dim. 2') +
  guides(color=F) +
  facet_wrap(~type, nrow=3, strip.position = 'right')
png('split-umap.png', height=3, width=2, res=300, units='in')
p
dev.off()
