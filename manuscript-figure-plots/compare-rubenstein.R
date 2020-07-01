#!/usr/bin/env Rscript
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)

args <- commandArgs(T)
CLUSTER_ASSIGNMENTS <- args[1]
RUBENSTEIN <- args[2]
OUT <- args[3]
COUNTS <-args[4:length(args)]

# get our Type I vs type II RNA-seq TPM
load_counts <- function(f) {
  tmp <- read.table(f, head=F, sep='\t', col.names=c('library', 'barcode', 'gene', 'count'), colClasses = c('character', 'character', 'character', 'numeric'))
  return(tmp)
}

counts <- bind_rows(lapply(COUNTS, load_counts))
clusters <- read.table(CLUSTER_ASSIGNMENTS, head=F, sep='\t', col.names = c('library', 'barcode', 'cluster'), colClasses = c('character', 'character', 'numeric'))
clusters <- left_join(clusters, counts)
clusters <- clusters[clusters$cluster %in% c(0, 1),] %>%
  dplyr::filter(!is.na(gene))

# Do type I / type II

fiber_type_counts <- clusters %>%
  dplyr::group_by(cluster, gene) %>%
  dplyr::summarize(count=sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(cluster)%>%
  dplyr::mutate(tpm=1e6*count/sum(count)) %>%
  dplyr::ungroup()

our_fold_changes <- fiber_type_counts %>%
  dplyr::mutate(cluster=ifelse(cluster==0, 'type_II', 'type_I')) %>%
  dplyr::select(cluster, gene, tpm) %>%
  tidyr::spread(key=cluster, value=tpm) %>%
  dplyr::mutate(our_log2FC=log2(type_II / type_I))

rubenstein <- read.table(RUBENSTEIN, head=T, sep='\t')
rubenstein$log2FC[rubenstein$fiber_type=='Type I'] <- -1*rubenstein$log2FC[rubenstein$fiber_type=='Type I']

both <- left_join(rubenstein %>% dplyr::select(gene, log2FC, padj), our_fold_changes %>% dplyr::select(gene, our_log2FC))
both$significant <- ifelse(both$padj<=0.05, 'Sign.', 'N.S.')
p <- ggplot(both[both$significant=='Sign.',]) +
  #geom_point(aes(x=log2FC, y=our_log2FC, color=significant)) +
  geom_point(aes(x=log2FC, y=our_log2FC)) +
  geom_text_repel(aes(x=log2FC, y=our_log2FC, label=gene)) +
  #scale_color_manual(values = c('Sign.'='red', 'N.S.'='black')) +
  #guides(color=guide_legend(title='Sign. in\nRubenstein')) +
  theme_bw() +
  xlab('Rubenstein et. al log2(Type IIa / Type I fibers)') +
  ylab('Our log2(Type II / Type I fibers)') +
  geom_hline(yintercept = 0, linetype='dashed', color='grey') +
  geom_vline(xintercept = 0, linetype='dashed', color='grey') +
  geom_abline(slope = 1, intercept = 0, linetype='dashed', color='black')
pdf('rubenstein-vs-our-fiber-type-lfcs.pdf', height=6, width=6)
p
dev.off()


