#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option(c("--nuclei"), action = "store", type = "character", default = "", help = "[Required] List of nuclei to keep (library, barcode)"),
  make_option(c("--atac-metrics"), dest='atac_metrics', action = "store", type = "character", default = "", help = "[Required] Ataqv metrics (library-barcode, metric, value)"),
  make_option(c("--rna-metrics"), dest='rna_metrics', action = "store", type = "character", default = "", help = "[Required] Read counts after filtering for RNA nuclei (library, barcode, count)"),
  make_option(c("--library-labels"), dest='library_labels', action = "store", type = "character", default = "", help = "[Required] path to library-labels.txt")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts <- parse_args(option_parser)

library(viridis)
library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)


nuclei <- read.table(opts$nuclei, head = F, as.is = T, sep = '\t', col.names = c('library', 'barcode'), colClasses = c('character', 'character'))
nuclei <- with(nuclei, paste(library, barcode, sep='-'))

atac_metrics <- read.table(opts$atac_metrics, head = F, as.is = T, sep = '\t', col.names = c('nucleus', 'metric', 'value')) %>%
  dplyr::filter(metric=='hqaa') %>%
  dplyr::select(nucleus, value) %>%
  dplyr::rename(hqaa=value) %>%
  dplyr::mutate(modality='ATAC',
                hqaa=as.numeric(hqaa)/2)

rna_metrics <- read.table(opts$rna_metrics, head = F, as.is = T, sep = '\t', col.names = c('library', 'barcode', 'count')) %>%
  dplyr::mutate(nucleus=glue('{library}-{barcode}')) %>%
  dplyr::rename(hqaa=count) %>%
  dplyr::select(nucleus, hqaa) %>%
  dplyr::mutate(modality='RNA')

all <- bind_rows(atac_metrics, rna_metrics) %>%
  dplyr::filter(nucleus %in% nuclei)
all$library <- gsub('(.*)-(.*)', '\\1', all$nucleus)



medians <- all %>%
  dplyr::group_by(library) %>%
  dplyr::summarize(med=median(hqaa),
                   avg=mean(hqaa))



# make plots for FANS analysis
fans <- all[grep('133', all$library),]
labels <- unique(fans[,c('library', 'modality')])
labels$fans_status <- ifelse(labels$library %in% c('133151-hg19', '133153-hg19', '133155-hg19', '133157-hg19'), 'no FANS', 'FANS')
labels$replicate <- ifelse(labels$library %in% c('133151-hg19', '133152-hg19', '133155-hg19', '133156-hg19'), 'rep. 1', 'rep. 2')
labels$label <- paste(labels$fans_status, labels$replicate, sep=', ')
labels <- labels[order(labels$fans_status, labels$replicate),]
labels$label <- factor(labels$label, levels=rev(unique(labels$label)), ordered = T)
fans <- left_join(fans, labels)

fans_atac <- fans[fans$modality=='ATAC',]
fans_rna <- fans[fans$modality=='RNA',]
p <- ggplot(fans_atac) +
  geom_jitter(aes(x = label, y = hqaa, color=modality), alpha = 0.3, stroke=0) +
  geom_errorbar(aes(x = label, ymin = med, ymax=med), data = left_join(fans_atac, medians) %>%
                  dplyr::select(modality, fans_status, replicate, med, avg, label) %>% unique()) +
  theme_bw() +
  scale_y_log10() +
  ylab('# fragments') +
  xlab('') +
  scale_color_viridis_d() +
  theme(axis.text.x = element_text(angle = 45, vjust=0.85, hjust=0.75)) +
  coord_flip() +
  guides(color=F)
png('fans-atac-nuclei-read-counts.png', height=3, width=4, units='in', res=300)
p
dev.off()

p <- ggplot(fans_rna) +
  geom_jitter(aes(x = label, y = hqaa, color=modality), alpha = 0.3, stroke=0) +
  geom_errorbar(aes(x = label, ymin = med, ymax=med), data = left_join(fans_rna, medians) %>%
                  dplyr::select(modality, fans_status, replicate, med, avg, label) %>% unique()) +
  theme_bw() +
  scale_y_log10() +
  ylab('# fragments') +
  xlab('') +
  scale_color_viridis_d(begin=1) +
  theme(axis.text.x = element_text(angle = 45, vjust=0.85, hjust=0.75)) +
  coord_flip() +
  guides(color=F)
png('fans-rna-nuclei-read-counts.pdf', height=5, width=7, units='in', res=300)
p
dev.off()

p <- fans %>%
  dplyr::group_by(label, modality) %>%
  dplyr::summarize(number_nuclei=n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(modality=='ATAC') %>%
  ggplot() +
  geom_bar(aes(x = label, y = number_nuclei, fill=modality), stat = 'identity') +
  geom_text(aes(x = label, y = number_nuclei+150, label=number_nuclei)) +
  theme_bw() +
  #ylim(c(0, max(counts$number_nuclei * 1.2))) +
  scale_fill_viridis_d() +
  ylab('# nuclei') +
  xlab('') +
  theme(axis.text.x = element_text(angle = 45, vjust=0.85, hjust=0.75)) +
  coord_flip() +
  guides(fill=F)
png('fans-atac-nuclei-per-library.png', height=5, width=5, units='in', res=300)
p
dev.off()

p <- fans %>%
  dplyr::group_by(label, modality) %>%
  dplyr::summarize(number_nuclei=n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(modality=='RNA') %>%
  ggplot() +
  geom_bar(aes(x = label, y = number_nuclei, fill=modality), stat = 'identity') +
  geom_text(aes(x = label, y = number_nuclei+500, label=number_nuclei)) +
  theme_bw() +
  #ylim(c(0, max(counts$number_nuclei * 1.2))) +
  scale_fill_viridis_d() +
  ylab('# nuclei') +
  xlab('') +
  theme(axis.text.x = element_text(angle = 45, vjust=0.85, hjust=0.75)) +
  coord_flip() +
  guides(fill=F)
png('fans-rna-nuclei-per-library.png', height=5, width=5, units='in', res=300)
p
dev.off()

# make plots for 20k vs 40k experiment
loading <- all[grep('63_', all$library),]
labels <- unique(loading[,c('library', 'modality')])
labels$label <- paste0(gsub('63_(\\d+).*', '\\1', labels$library), 'k')
loading <- left_join(loading, labels)

loading_atac <- loading[loading$modality=='ATAC',]
loading_rna <- loading[loading$modality=='RNA',]
p <- ggplot(loading_atac) +
  geom_jitter(aes(x = label, y = hqaa, color=modality), alpha = 0.3, stroke=0) +
  geom_errorbar(aes(x = label, ymin = med, ymax=med), data = left_join(loading_atac, medians) %>%
                  dplyr::select(modality, med, avg, label) %>% unique()) +
  theme_bw() +
  scale_y_log10() +
  ylab('# fragments') +
  xlab('') +
  scale_color_viridis_d() +
  theme(axis.text.x = element_text(angle = 45, vjust=0.85, hjust=0.75)) +
  coord_flip() +
  guides(color=F)
png('loading-atac-nuclei-read-counts.png', height=3, width=4, units='in', res=300)
p
dev.off()

p <- ggplot(loading_rna) +
  geom_jitter(aes(x = label, y = hqaa, color=modality), alpha = 0.3, stroke=0) +
  geom_errorbar(aes(x = label, ymin = med, ymax=med), data = left_join(loading_rna, medians) %>%
                  dplyr::select(modality, med, avg, label) %>% unique()) +
  theme_bw() +
  scale_y_log10() +
  ylab('# fragments') +
  xlab('') +
  scale_color_viridis_d(begin=1) +
  theme(axis.text.x = element_text(angle = 45, vjust=0.85, hjust=0.75)) +
  coord_flip() +
  guides(color=F)
png('loading-rna-nuclei-read-counts.pdf', height=5, width=7, units='in', res=300)
p
dev.off()

p <- loading %>%
  dplyr::group_by(label, modality) %>%
  dplyr::summarize(number_nuclei=n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(modality=='ATAC') %>%
  ggplot() +
  geom_bar(aes(x = label, y = number_nuclei, fill=modality), stat = 'identity') +
  geom_text(aes(x = label, y = number_nuclei+150, label=number_nuclei)) +
  theme_bw() +
  #ylim(c(0, max(counts$number_nuclei * 1.2))) +
  scale_fill_viridis_d() +
  ylab('# nuclei') +
  xlab('') +
  theme(axis.text.x = element_text(angle = 45, vjust=0.85, hjust=0.75)) +
  coord_flip() +
  guides(fill=F)
png('loading-atac-nuclei-per-library.png', height=5, width=5, units='in', res=300)
p
dev.off()

p <- loading %>%
  dplyr::group_by(label, modality) %>%
  dplyr::summarize(number_nuclei=n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(modality=='RNA') %>%
  ggplot() +
  geom_bar(aes(x = label, y = number_nuclei, fill=modality), stat = 'identity') +
  geom_text(aes(x = label, y = number_nuclei+500, label=number_nuclei)) +
  theme_bw() +
  #ylim(c(0, max(counts$number_nuclei * 1.2))) +
  scale_fill_viridis_d() +
  ylab('# nuclei') +
  xlab('') +
  theme(axis.text.x = element_text(angle = 45, vjust=0.85, hjust=0.75)) +
  coord_flip() +
  guides(fill=F)
png('loading-rna-nuclei-per-library.png', height=5, width=5, units='in', res=300)
p
dev.off()


# get summary stats for all the quality data...
library_labels <- read.table(opts$library_labels, head = T, as.is = T, sep='\t')
summary_stats <- all %>%
  dplyr::group_by(library) %>%
  dplyr::summarize(number_nuclei=n(),
                   mean_fragments=mean(hqaa),
                   median_fragments=median(hqaa)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!library %in% c('133152-hg19', '133154-hg19')) %>%
  left_join(library_labels)
summary_stats <- summary_stats[,c('name', 'species', 'modality', 'replicate', 'samples', 'experiment', 'number_nuclei', 'mean_fragments', 'median_fragments')]
write.table(summary_stats, file = 'summary_stats.tsv', append = F, quote = F, sep = '\t', row.names = F, col.names = T)
