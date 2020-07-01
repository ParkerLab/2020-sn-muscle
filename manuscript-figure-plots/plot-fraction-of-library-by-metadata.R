#!/usr/bin/env Rscript
# plot the cell type representation per cluster
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(glue)

args <- commandArgs(T)
CLUSTER_FILE <- args[1]
CLUSTER_NAMES_FILE <- args[2]
NUCLEI_WITH_INDIVIDUAL_ASSIGNMENTS <- args[3]
LIBRARY_LABELS <- args[4]

modalities <- read.table(LIBRARY_LABELS, sep='\t', head=T, as.is=T)[,c('library', 'modality')]
individuals <- read.table(NUCLEI_WITH_INDIVIDUAL_ASSIGNMENTS, head=F, sep='\t', as.is=T, col.names=c('library', 'barcode', 'ind'))
individuals$ind[individuals$ind=='rat1'] <- 'rat'

clusters <- read.table(CLUSTER_FILE, head = F, sep = '\t', col.names = c('library', 'barcode', 'cluster'), colClasses = c('character', 'character', 'numeric'))
cluster_names <- read.table(CLUSTER_NAMES_FILE, head = T, sep = '\t') %>% rename(cluster=old_name, cluster_name=new_name)
clusters <- left_join(clusters, individuals)
clusters$label <- with(clusters, paste0(library, ', ', ind))
clusters <- left_join(clusters, cluster_names) %>%
  dplyr::group_by(library, label, cluster, cluster_name, ind) %>%
  dplyr::summarize(count=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(library, label, ind) %>%
  dplyr::mutate(fraction_of_library=count/sum(count)) %>%
  ungroup() %>%
  left_join(modalities)


# first, look at consistency within individuals across libraries...
# generally
means <- clusters %>%
  dplyr::group_by(ind, modality, cluster_name) %>%
  dplyr::summarize(mean_fraction=mean(fraction_of_library))
  
p <- ggplot(clusters) +
  geom_bar(aes(x = library, y = fraction_of_library, fill = modality), position='dodge', stat='identity') +
  geom_text(aes(x = library, y = fraction_of_library*1.1, label = count), position='dodge') +
  geom_hline(aes(yintercept = mean_fraction, color=modality), data=means, linetype='dashed') +
  ylab('Fraction of nuclei') +
  facet_grid(cluster_name~ind, scales='free', space = 'free_x') +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  xlab('Library')
pdf('fraction-of-library-within-individual.barplot.pdf', width=15, height=15)
p
dev.off()

p <- ggplot(clusters) +
  geom_jitter(aes(x = ind, y = fraction_of_library, color = modality), height=0, width=0.2) +
  geom_errorbar(aes(x = ind, ymin = mean_fraction, ymax = mean_fraction, color=modality), data=means, linetype='dashed') +
  ylab('Fraction of nuclei') +
  facet_wrap(~cluster_name, scales='free_y') +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
  scale_color_viridis_d() +
  xlab('Biological sample')
pdf('fraction-of-library-within-individual.jitter.pdf', width=8.5, height=6)
p
dev.off()

# chisq test on KSM1 to compare modalities?



# compare to Seurat
