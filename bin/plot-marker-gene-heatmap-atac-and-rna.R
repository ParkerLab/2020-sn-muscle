#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(glue)

args <- commandArgs(T)
COUNTS <- args[1]
CLUSTERS <- args[2]

counts <- read.table(COUNTS, head = F, sep = '\t', col.names = c('cluster', 'modality', 'feature', 'count'), colClasses = c('character', 'character', 'character', 'numeric'))
counts$species <- ifelse(gsub('.*-', '', counts$cluster)=='hg19', 'human', 'rat')
counts$cluster <- as.numeric(gsub('-.*', '', counts$cluster))
counts <- counts %>%
  dplyr::group_by(species, modality, cluster) %>%
  dplyr::mutate(tpm=1e6*count/sum(count)) %>%
  dplyr::ungroup()
counts <- counts %>%
  dplyr::group_by(species, modality, feature) %>%
  dplyr::mutate(relative_tpm=tpm/max(tpm)) %>%
  dplyr::ungroup()

both <- counts
FEATURE_ORDERS <- c('MYH1', 'MYH7', 'PDGFRA', 'VWF', 'MYH11', 'CD163', 'PAX7')
both <- both[both$feature %in% FEATURE_ORDERS,]
both$feature <- factor(both$feature, levels=FEATURE_ORDERS, ordered=T)
both$type <- with(both, paste(species, modality))
both$cluster <- factor(both$cluster, levels=rev(sort(unique(both$cluster))), ordered=T)

tmp <- both %>% 
  complete(feature, cluster, type, fill=list('relative_tpm'=0))
tmp$type <- factor(tmp$type, levels = c('human RNA', 'human ATAC', 'rat ATAC'), ordered = T)
p <- ggplot(tmp) +
  geom_tile(aes(x = feature, y = cluster, fill=relative_tpm)) +
  theme_bw() +
  xlab('Gene') +
  ylab('Cluster') +
  facet_wrap(~type, nrow=1) +
  scale_fill_gradient(low='white', high='red') +
  guides(fill=guide_colorbar(title='Rel. expression\nor accessibility')) +
  theme(legend.position = 'top', axis.text.x = element_text(angle=45, vjust=1, hjust=1))# +
pdf('marker-genes-rat-and-human.pdf', height=3, width=4.5)
p
dev.off()





both <- both[both$species=='human',]

both <- both %>% 
  complete(feature, cluster, type, fill=list('relative_tpm'=0))
p <- ggplot(both) +
  geom_tile(aes(x = feature, y = cluster, fill=relative_tpm)) +
  theme_bw() +
  xlab('Gene') +
  ylab('Cluster') +
  facet_wrap(~type, nrow=1) +
  scale_fill_gradient(low='white', high='red') +
  guides(fill=guide_colorbar(title='Rel. expression\nor accessibility')) +
  theme(legend.position = 'top', axis.text.x = element_text(angle=45, vjust=1, hjust=1))# +
pdf('marker-genes.pdf', height=3, width=3.5)
p
dev.off()
