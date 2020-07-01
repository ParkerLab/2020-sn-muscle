#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(GGally)

COUNTS <- commandArgs(T)[1]
counts <- read.table(COUNTS, head = T, as.is = T, sep = '\t', col.names = c('library', 'feature', 'count'))
counts <- counts[counts$feature!='nan',]

library_labels <- read.table('library-labels.txt', head = T, as.is = T, sep = '\t')
library_labels <- library_labels[library_labels$modality=='RNA',]
fans_info <- library_labels[!is.na(library_labels$experiment) & library_labels$experiment=='FANS vs. no FANS',]
#fans_info$fans_status <- ifelse(grepl('no FANS', fans_info$name), 'Crude', 'FANS')
fans_info$fans_status <- ifelse(grepl('no FANS', fans_info$name), 'no FANS', 'FANS')
fans_info$label <- paste(fans_info$fans_status, '\n', fans_info$replicate)
loading_info <- library_labels[!is.na(library_labels$experiment) & library_labels$experiment=='20k vs. 40k',]
loading_info$loading_concentration <- ifelse(grepl('20', loading_info$name), '20k', '40k')
loading_info$label <- loading_info$loading_concentration
all_info <- bind_rows(fans_info[,c('library', 'label')], loading_info[,c('library', 'label')])

# plot correlation for FANS
fans <- counts[counts$library %in% fans_info$library,] %>%
  left_join(fans_info[,c('library', 'label')]) %>%
  dplyr::select(-library) %>%
  dplyr::group_by(label) %>%
  dplyr::mutate(tpm=1e6*count/sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(log_tpm=log10(tpm+0.01)) %>%
  dplyr::select(label, feature, log_tpm) %>%
  tidyr::spread(key=label, value=log_tpm, fill=log10(0.01))


p <- ggpairs(fans, columns = 2:ncol(fans), lower = list(continuous = wrap("points", alpha = 0.1, stroke=0))) +
  theme_bw() +
  xlab('log10(CPM + 0.01)') +
  ylab('log10(CPM + 0.01)')
png('fans-rna-correlation.png', width = 3.3, height = 3.3, units='in', res=300)
p
dev.off()


# plot correlation for 20k vs 40k
loading <- counts[counts$library %in% loading_info$library,] %>%
  left_join(loading_info[,c('library', 'label')]) %>%
  dplyr::select(-library) %>%
  dplyr::group_by(label) %>%
  dplyr::mutate(tpm=1e6*count/sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(log_tpm=log10(tpm+0.01)) %>%
  dplyr::select(label, feature, log_tpm) %>%
  tidyr::spread(key=label, value=log_tpm, fill=log10(0.01))


p <- ggpairs(loading, columns = 2:ncol(loading), lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw() +
  xlab('log10(CPM + 0.01)') +
  ylab('log10(CPM + 0.01)')
png('loading-rna-correlation.png', width = 2.3, height = 2.3, units='in', res=300)
p
dev.off()


# correlation for all RNA-seq
all <- counts[counts$library %in% all_info$library,] %>%
  left_join(all_info[,c('library', 'label')]) %>%
  dplyr::select(-library) %>%
  dplyr::group_by(label) %>%
  dplyr::mutate(tpm=1e6*count/sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(log_tpm=log10(tpm+0.01)) %>%
  dplyr::select(label, feature, log_tpm) %>%
  tidyr::spread(key=label, value=log_tpm, fill=log10(0.01))


p <- ggpairs(all, columns = 2:ncol(all), lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw() +
  xlab('log10(CPM + 0.01)') +
  ylab('log10(CPM + 0.01)')
png('all-rna-correlation.png', width = 7, height = 7, units='in', res=300)
p
dev.off()
