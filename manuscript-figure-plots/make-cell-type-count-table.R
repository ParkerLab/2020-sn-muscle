#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)

args <- commandArgs(T)

CLUSTER_NAMES <- args[1]
CLUSTER_ASSIGNMENTS <- args[2]
LIBRARY_LABELS <- args[3]

library_to_modality_and_species <- read.table(LIBRARY_LABELS, head=T, as.is=T, sep='\t') %>%
  dplyr::select(library, species, modality)

clusters <- read.table(CLUSTER_ASSIGNMENTS, head=F, sep='\t', col.names = c('library', 'barcode', 'cluster'), colClasses = c('character'))
cluster_names <- read.table(CLUSTER_NAMES, head=T, sep='\t', colClasses = c('character')) %>% dplyr::rename(cluster=old_name)
clusters <- left_join(clusters, cluster_names)
clusters <- left_join(clusters, library_to_modality_and_species)
clusters$label <- paste(clusters$species, clusters$modality, sep='_')
tbl <- clusters %>%
  dplyr::group_by(label, new_name) %>%
  dplyr::summarize(count=n()) %>%
  dplyr::ungroup() %>%
  tidyr::spread(key=label, value=count) %>%
  dplyr::rename(cell_type=new_name)
tbl <- tbl[rev(order(tbl$Human_RNA)),]

write.table(tbl, file = 'cell_type_counts.tsv', append = F, quote = F, sep = '\t', row.names = F, col.names = T)
