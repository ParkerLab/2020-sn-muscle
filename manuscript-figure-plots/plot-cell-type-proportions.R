#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(tidyr)

args <- commandArgs(T)

CLUSTER_NAMES <- args[1]
CLUSTER_ASSIGNMENTS <- args[2]

clusters <- read.table(CLUSTER_ASSIGNMENTS, head=F, sep='\t', col.names = c('library', 'barcode', 'cluster'), colClasses = c('character'))
cluster_names <- read.table(CLUSTER_NAMES, head=T, sep='\t', colClasses = c('character')) %>% dplyr::rename(cluster=old_name)
clusters <- left_join(clusters, cluster_names)
clusters <- clusters %>% dplyr::group_by(new_name, cluster) %>%
  dplyr::summarize(count=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(frac=count/sum(count))
clusters <- clusters[order(clusters$cluster),]
clusters$new_name <- factor(clusters$new_name, levels=rev(clusters$new_name), ordered=T)
labels <- clusters
labels$pos <- labels$frac
labels$pos[100*labels$pos<5] <- labels$pos[100*labels$pos<5] + (20/100)
labels$pos[100*labels$pos>5] <- labels$pos[100*labels$pos>5] - (10/100)
labels$label <- paste0(round(100*labels$frac, 1), '%')
labels$colors <- ifelse(100*labels$frac<5, 'black', 'white')

p <- ggplot(clusters) +
  geom_bar(aes(x=new_name, y=frac, fill=new_name), stat='identity') +
  geom_text(aes(x=new_name, y=pos, label=label, color=colors), data=labels) +
  theme_bw() +
  coord_flip() +
  scale_fill_viridis_d(direction = -1) +
  scale_color_manual(values=rev(unique(labels$colors))) +
  guides(fill=F, color=F) +
  ylab('Fraction of nuclei') +
  xlab('')
pdf('cell-type-proportions.pdf', height=1.7, width=3)
p
dev.off()
