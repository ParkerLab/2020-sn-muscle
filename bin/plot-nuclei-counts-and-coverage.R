#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option(c("--nuclei"), action = "store", type = "character", default = "", help = "[Required] List of nuclei to keep (library, barcode)"),
  make_option(c("--atac-metrics"), dest='atac_metrics', action = "store", type = "character", default = "", help = "[Required] Ataqv metrics (library-barcode, metric, value)"),
  make_option(c("--rna-metrics"), dest='rna_metrics', action = "store", type = "character", default = "", help = "[Required] Read counts after filtering for RNA nuclei (library, barcode, count)")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts <- parse_args(option_parser)

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
                hqaa=as.numeric(hqaa))

rna_metrics <- read.table(opts$rna_metrics, head = F, as.is = T, sep = '\t', col.names = c('library', 'barcode', 'count')) %>%
  dplyr::mutate(nucleus=glue('{library}-{barcode}')) %>%
  dplyr::rename(hqaa=count) %>%
  dplyr::select(nucleus, hqaa) %>%
  dplyr::mutate(modality='RNA')

all <- bind_rows(atac_metrics, rna_metrics) %>%
  dplyr::filter(nucleus %in% nuclei)
all$library <- gsub('(.*)-(.*)', '\\1', all$nucleus)

# plot the read counts for each nucleus

medians <- all %>%
  dplyr::group_by(library) %>%
  dplyr::summarize(med=median(hqaa))

p <- ggplot(all) +
  geom_jitter(aes(x = library, y = hqaa, color=modality), alpha = 0.3, stroke=0) +
  geom_errorbar(aes(x = library, ymin = med, ymax=med), data = medians) +
  theme_bw() +
  scale_color_viridis_d() +
  scale_y_log10() +
  ylab('# reads') +
  xlab('') +
  theme(axis.text.x = element_text(angle = 45, vjust=0.85, hjust=0.75)) +
  coord_flip()
pdf('nuclei-read-counts.pdf', height=5, width=7)
p
dev.off()


# plot the number of nuclei for each library

counts <- all %>%
  dplyr::group_by(library, modality) %>%
  dplyr::summarize(number_nuclei=n()) %>%
  dplyr::mutate(label=glue('{number_nuclei}'))

p <- ggplot(counts) +
  geom_bar(aes(x = library, y = number_nuclei, fill=modality), stat = 'identity') +
  geom_text(aes(x = library, y = number_nuclei+1000, label=label)) +
  theme_bw() +
  ylim(c(0, max(counts$number_nuclei * 1.2))) +
  scale_fill_viridis_d() +
  ylab('# nuclei') +
  xlab('') +
  theme(axis.text.x = element_text(angle = 45, vjust=0.85, hjust=0.75)) +
  coord_flip()
pdf('nuclei-per-library.pdf', height=5, width=5)
p
dev.off()

# plot number nuclei per modality
counts <- all %>%
  dplyr::group_by(modality) %>%
  dplyr::summarize(number_nuclei=n()) %>%
  dplyr::mutate(label=glue('{number_nuclei}'))

p <- ggplot(counts) +
  geom_bar(aes(x = modality, y = number_nuclei, fill=modality), stat = 'identity') +
  geom_text(aes(x = modality, y = number_nuclei+1000, label=label)) +
  theme_bw() +
  ylim(c(0, max(counts$number_nuclei * 1.2))) +
  scale_fill_viridis_d() +
  ylab('# nuclei') +
  xlab('') +
  theme(axis.text.x = element_text(angle = 45, vjust=0.85, hjust=0.75)) +
  coord_flip() +
  guides(fill=F)
pdf('nuclei-per-modality.pdf', height=2, width=4)
p
dev.off()

