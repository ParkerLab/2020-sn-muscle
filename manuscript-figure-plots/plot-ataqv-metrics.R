#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)

library_labels <- read.table('library-labels.txt', head = T, as.is = T, sep = '\t')
#library_labels$fans_status <- ifelse(grepl('no FANS', library_labels$fans_status), 'Crude', 'FANS')
library_labels$fans_status <- ifelse(grepl('no FANS', library_labels$fans_status), 'no FANS', 'FANS')
library_labels$loading_concentration <- NA
library_labels$loading_concentration[grepl('20', library_labels$name)] <- '20k'
library_labels$loading_concentration[grepl('40', library_labels$name)] <- '40k'

metrics <- read.table('metrics.txt', head = F, as.is = T, sep = '\t', col.names=c('library', 'metric', 'value'))
metrics$metric[metrics$metric=='total_peaks'] <- 'Num. peaks'
metrics$metric[metrics$metric=='hqaa_overlapping_peaks_percent'] <- 'Perc. reads overlapping peaks'
metrics <- left_join(metrics, library_labels)

# plots for FANS experiment first
bulk_info <- data.frame(library=c('320-NM-1', '320-NM-2'), label=c('Bulk, rep. 1', 'Bulk, rep. 2'))
p <- metrics %>%
  dplyr::filter(!is.na(experiment) & experiment=="FANS vs. no FANS") %>%
  dplyr::mutate(label=glue('{fans_status}, {replicate}')) %>% 
  ggplot() +
  geom_point(aes(x = fans_status, y=value), alpha=0.5) +
  facet_wrap(~metric, scales = 'free') +
  theme_bw() +
  xlab('')
pdf('fans-aggregate-ataqv-metrics.pdf', height=2, width=6)
p
dev.off()

#COLORS <- c('Crude, rep. 1'='#8c510a', 'Crude, rep. 2'='#d8b365', 'FANS, rep. 1'='#5ab4ac', 'FANS, rep. 2'='#01665e', 'Bulk, rep. 1'='black', 'Bulk, rep. 2'='grey')
COLORS <- c('no FANS, rep. 1'='#8c510a', 'no FANS, rep. 2'='#d8b365', 'FANS, rep. 1'='#5ab4ac', 'FANS, rep. 2'='#01665e', 'Bulk, rep. 1'='black', 'Bulk, rep. 2'='grey')
flds <- read.table('fragment-lengths.txt', head=F, as.is = T, sep = '\t', col.names = c('library', 'fragment_length', 'count')) %>%
  group_by(library) %>%
  dplyr::mutate(proportion=count/sum(count)) %>%
  dplyr::ungroup() %>%
  left_join(library_labels) %>%
  dplyr::filter(!is.na(experiment) & experiment=="FANS vs. no FANS") %>%
  dplyr::mutate(label=glue('{fans_status}, {replicate}')) %>%
  dplyr::select(fragment_length, proportion, label)
flds2 <- read.table('fragment-lengths.txt', head=F, as.is = T, sep = '\t', col.names = c('library', 'fragment_length', 'count')) %>%
  group_by(library) %>%
  dplyr::mutate(proportion=count/sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(library %in% bulk_info$library) %>%
  left_join(bulk_info) %>%
  dplyr::select(fragment_length, proportion, label)
flds <- bind_rows(flds, flds2)
p <- ggplot(flds) +
  geom_line(aes(x = fragment_length, y=proportion, color=label)) +
  theme_bw() +
  #scale_color_viridis_d(option='C') +
  #scale_color_brewer(palette='BrBG') +
  scale_color_manual(values=COLORS) +
  xlab('Fragment length (bp)') +
  ylab('Fraction of fragments') +
  xlim(c(0, 500)) +
  guides(color=guide_legend(title='', ncol = 1)) +
  #guides(color=guide_legend(title=''))
  #guides(color=F) +
  #theme(plot.margin = unit(c(0, 0, 0, 0), 'cm')) +
  theme(legend.position = 'right', legend.margin = margin(0, 0, 0, 0), legend.text = element_text(size=10), legend.box.margin = margin(-1, 0, -10, 0), plot.margin = unit(c(0, 0, 0, 0), 'cm'))
pdf('fans-aggregate-fld.pdf', height=2, width=4)
p
dev.off()


tss <- read.table('tss-coverage.txt', head=F, as.is = T, sep = '\t', col.names = c('library', 'position', 'value')) %>%
  left_join(library_labels) %>%
  dplyr::filter(!is.na(experiment) & experiment=="FANS vs. no FANS") %>%
  dplyr::mutate(label=glue('{fans_status}, {replicate}')) %>%
  dplyr::select(position, value, label)
tss2 <- read.table('tss-coverage.txt', head=F, as.is = T, sep = '\t', col.names = c('library', 'position', 'value')) %>%
  left_join(library_labels) %>%
  dplyr::filter(library %in% bulk_info$library) %>%
  left_join(bulk_info) %>%
  dplyr::select(position, value, label)
tss <- bind_rows(tss, tss2)
p <- ggplot(tss) +
  geom_line(aes(x = position, y=value, color=label)) +
  theme_bw() +
  #scale_color_viridis_d(option='C') +
  #scale_color_brewer(palette='BrBG') +
  scale_color_manual(values=COLORS) +
  xlab('Pos. relative to TSS') +
  ylab('Normalized coverage') +
  guides(color=guide_legend(title='', ncol = 1)) +
  #guides(color=guide_legend(title='', nrow = 2)) +
  #theme(legend.position = 'top', legend.margin = margin(0, 1, 0, 0), legend.box.margin = margin(-1, 1, -10, 0), plot.margin = unit(c(0, 1, 0, 0), 'cm'))
  theme(legend.position = 'right', legend.margin = margin(0, 0, 0, 0), legend.text = element_text(size=10), legend.box.margin = margin(-1, 0, -10, 0), plot.margin = unit(c(0, 0, 0, 0), 'cm'))
pdf('fans-aggregate-tss.pdf', height=2, width=4)
p
dev.off()


# plots for loading experiments
bulk_info <- data.frame(library=c('320-NM-1', '320-NM-2', '320-NM-3', '320-NM-4'), label=c('HSM1 bulk, rep. 1', 'HSM1 bulk, rep. 2', 'HSM2 bulk, rep. 1', 'HSM2 bulk, rep. 2'))
bulk_info$label <- gsub(', ', '\n', bulk_info$label)
COLORS <- c('40k'='#8c510a', '20k'='#d8b365', 'HSM2 bulk, rep. 1'='#5ab4ac', 'HSM2 bulk, rep. 2'='#01665e', 'HSM1 bulk, rep. 1'='black', 'HSM1 bulk, rep. 2'='grey')
names(COLORS) <- gsub(', ', '\n', names(COLORS))
#COLORS <- c('20k'='#762a83', '40k'='#1b7837')
flds <- read.table('fragment-lengths.txt', head=F, as.is = T, sep = '\t', col.names = c('library', 'fragment_length', 'count')) %>%
  group_by(library) %>%
  dplyr::mutate(proportion=count/sum(count)) %>%
  dplyr::ungroup() %>%
  left_join(library_labels) %>%
  dplyr::filter(!is.na(experiment) & experiment=="20k vs. 40k") %>%
  dplyr::mutate(label=glue('{loading_concentration}')) %>%
  dplyr::select(fragment_length, proportion, label)
flds2 <- read.table('fragment-lengths.txt', head=F, as.is = T, sep = '\t', col.names = c('library', 'fragment_length', 'count')) %>%
  group_by(library) %>%
  dplyr::mutate(proportion=count/sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(library %in% bulk_info$library) %>%
  left_join(bulk_info) %>%
  dplyr::select(fragment_length, proportion, label)
flds <- bind_rows(flds, flds2)
p <- ggplot(flds) +
  geom_line(aes(x = fragment_length, y=proportion, color=label)) +
  theme_bw() +
  #scale_color_viridis_d(option='C') +
  #scale_color_brewer(palette='PRGn') +
  scale_color_manual(values=COLORS) +
  xlab('Fragment length (bp)') +
  ylab('Fraction of fragments') +
  xlim(c(0, 500)) +
  #guides(color=guide_legend(title=''))
  #guides(color=F) +
  guides(color=guide_legend(title='', ncol = 1)) +
  #theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'))
  theme(legend.position = 'right', legend.margin = margin(0, 0, 0, 0), legend.text = element_text(size=11), legend.box.margin = margin(-1, 0, -10, 0), plot.margin = unit(c(0, 0, 0, 0), 'cm'))
pdf('loading-aggregate-fld.pdf', height=2, width=4)
p
dev.off()


tss <- read.table('tss-coverage.txt', head=F, as.is = T, sep = '\t', col.names = c('library', 'position', 'value')) %>%
  left_join(library_labels) %>%
  dplyr::filter(!is.na(experiment) & experiment=="20k vs. 40k") %>%
  dplyr::mutate(label=glue('{loading_concentration}')) %>%
  dplyr::select(position, value, label)
tss2 <- read.table('tss-coverage.txt', head=F, as.is = T, sep = '\t', col.names = c('library', 'position', 'value')) %>%
  left_join(library_labels) %>%
  dplyr::filter(library %in% bulk_info$library) %>%
  left_join(bulk_info) %>%
  dplyr::select(position, value, label)
tss <- bind_rows(tss, tss2)

p <- ggplot(tss) +
  geom_line(aes(x = position, y=value, color=label)) +
  theme_bw() +
  #scale_color_viridis_d(option='C') +
  #scale_color_brewer(palette='BrBG') +
  scale_color_manual(values=COLORS) +
  xlab('Pos. relative to TSS') +
  ylab('Normalized coverage') +
  #guides(color=guide_legend(title='')) +
  guides(color=guide_legend(title='', ncol = 1)) +
  #theme(legend.position = 'top', legend.margin = margin(0, 1, 0, 0), legend.box.margin = margin(-1, 1, -10, 0), plot.margin = unit(c(0, 1, 0, 0), 'cm'))
  theme(legend.position = 'right', legend.margin = margin(0, 0, 0, 0), legend.text = element_text(size=11), legend.box.margin = margin(-1, 0, -10, 0), plot.margin = unit(c(0, 0, 0, 0), 'cm'))
pdf('loading-aggregate-tss.pdf', height=2, width=4)
p
dev.off()
