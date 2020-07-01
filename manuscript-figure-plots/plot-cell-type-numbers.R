#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(tidyr)
library(glue)

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

clusters <- clusters %>% dplyr::group_by(new_name, cluster, label) %>%
  dplyr::summarize(count=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(frac=count/sum(count))
clusters <- clusters[order(clusters$cluster),]
clusters$new_name <- factor(clusters$new_name, levels=rev(unique(clusters$new_name)), ordered=T)
labels <- clusters
labels$pos <- labels$count*1.1
labels$colors <- ifelse(100*labels$frac<5, 'black', 'white')

for(LABEL in unique(labels$label)) {
  tmp <- labels[labels$label==LABEL,]
  tmp$pos <- tmp$count + max(tmp$count) * 0.1
  tmp$colors <- ifelse(tmp$count + max(tmp$count) * 0.1==max(tmp$pos), 'white', 'black')
  tmp$pos[tmp$pos==max(tmp$pos)] <- tmp$count[tmp$pos==max(tmp$pos)] * 0.8
  p <- ggplot(tmp) +
    geom_bar(aes(x=new_name, y=count, fill=new_name), stat='identity') +
    geom_text(aes(x=new_name, y=pos, label=count, color=colors), data=tmp) +
    theme_bw() +
    coord_flip() +
    scale_fill_viridis_d(direction = -1) +
    scale_color_manual(values=rev(unique(labels$colors))) +
    guides(fill=F, color=F) +
    ylab('# snATAC-seq nuclei') +
    xlab('') +
    scale_x_discrete(position='top') +
    theme(plot.margin =  margin(0, 0, 0, 0, "cm"))
  pdf(glue('cell-type-proportions-{LABEL}.pdf'), height=1.7, width=3.2)
  print(p)
  dev.off()
}
