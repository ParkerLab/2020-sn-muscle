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

# for testing
#DIM_FILE <- '/lab/work/porchard/sn-muscle-project/work/liger/results/liger-visualize/factorization_k_20___factorization_lambda_10___norm_knnk_20___norm_resolution_1___umap_neighbors_30.dim.txt'
#CLUSTER_FILE <- '/lab/work/porchard/sn-muscle-project/work/liger/results/liger-normalize/factorization_k_20___factorization_lambda_10___norm_knnk_20___norm_resolution_1.clusters.txt'

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


CLUSTERS <- unique(tmp$cluster)
facetted <- bind_rows(lapply(CLUSTERS, function(x){y <- tmp; y$color_cluster <- x; y$col <- ifelse(y$cluster==x, 'yes', 'no'); return(y)}))

colors <- c('yes'='red', 'no'='grey')
alphas <- c('yes'=1, 'no'=0.3)

p <- ggplot(facetted) +
  geom_point(aes(x = dim1, y = dim2, color = col, alpha=col), size = 1, stroke=0) +
  ylab('Dim. 2') +
  xlab('Dim. 1') +
  theme_bw() +
  scale_color_manual(values=colors) +
  scale_alpha_manual(values=alphas) +
  guides(color = guide_legend(title='Cluster')) +
  facet_wrap(~color_cluster)
png(OUT, height = 3*sqrt(length(CLUSTERS)), width = 3*sqrt(length(CLUSTERS)), res=300, units='in')
p
dev.off()
