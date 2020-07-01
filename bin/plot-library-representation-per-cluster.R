#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)

args <- commandArgs(T)
CLUSTER_FILE <- args[1]
PREFIX <- args[2]

# for testing
# CLUSTER_FILE <- '/lab/work/porchard/sn-muscle-project/work/downstream/results/liger/round-1/second-louvain-clusters.txt'

clusters <- read.table(CLUSTER_FILE, head = F, sep = '\t', col.names = c('library', 'barcode', 'cluster'), colClasses = c('character', 'character', 'numeric'))

clusters <- clusters %>%
	dplyr::group_by(cluster) %>%
	dplyr::mutate(size=n(),
		      label=glue('{cluster} (n = {size})'))
clusters <- clusters[rev(order(clusters$size)),]
clusters$label <- factor(clusters$label, levels=unique(clusters$label), ordered=T)

clusters <- clusters %>%
  dplyr::select(library, barcode, label) %>%
  dplyr::group_by(library, label) %>%
  dplyr::summarize(count=n())

p <- ggplot(clusters) +
  geom_point(aes(x = label, y = count, color = library)) +
  geom_line(aes(x = label, y = count, color = library, group=library)) +
  ylab('# nuclei') +
  xlab('Cluster') +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90)) +
  scale_color_viridis_d() +
  guides(color = guide_legend(title='Library'))
pdf(glue('{PREFIX}number-nuclei-per-cluster-per-library.pdf'), height = 5, width = 6)
p
dev.off()

tmp <- clusters %>%
  dplyr::group_by(library) %>%
  dplyr::mutate(fraction_of_library=count/sum(count))

p <- ggplot(tmp) +
  geom_point(aes(x = label, y = fraction_of_library, color = library)) +
  geom_line(aes(x = label, y = fraction_of_library, color = library, group=library)) +
  ylab('Fraction of library') +
  xlab('Cluster') +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90)) +
  scale_color_viridis_d() +
  guides(color = guide_legend(title='Library'))
pdf(glue('{PREFIX}fraction-library-per-cluster.pdf'), height = 5, width = 6)
p
dev.off()



tmp <- clusters %>%
  dplyr::group_by(label) %>%
  dplyr::mutate(fraction_of_cluster=count/sum(count))


p <- ggplot(tmp) +
  geom_point(aes(x = label, y = fraction_of_cluster, color = library)) +
  geom_line(aes(x = label, y = fraction_of_cluster, color = library, group=library)) +
  ylab('Fraction of cluster') +
  xlab('Cluster') +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90)) +
  scale_color_viridis_d() +
  guides(color = guide_legend(title='Library'))
pdf(glue('{PREFIX}fraction-cluster-per-library.pdf'), height = 5, width = 6)
p
dev.off()


