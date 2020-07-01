#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)
library(ggrepel)

args <- commandArgs(T)
CLUSTER_NAMES <- args[1]
ANNOTATED_LOCI <- args[2]

cluster_names <- read.table(CLUSTER_NAMES, head=T, sep='\t') %>%
  mutate(cluster=glue('cluster_{old_name}')) %>%
  dplyr::select(-old_name)
cluster_names$new_name <- as.character(cluster_names$new_name)
cluster_names$cluster <- as.character(cluster_names$cluster)

# load in the annotations
annotations <- read.table(ANNOTATED_LOCI, head=T, as.is=T, sep='\t')
colnames(annotations)[1:3] <- c('chrom', 'pos', 'locus')
annotations$locus <- as.character(annotations$locus)

# make a table:
df <- annotations %>%
  tidyr::gather(key=annotation, value=overlap, -chrom, -pos, -locus)
df$tissue <- NA
df$tissue[grep('adipose', df$annotation, ignore.case = T)] <- 'Adipose'
df$tissue[grep('beta', df$annotation, ignore.case = T)] <- 'Beta cell'
df$tissue[grep('liver', df$annotation, ignore.case = T)] <- 'Liver'
df$tissue[grep('islet', df$annotation, ignore.case = T)] <- 'Islets'
df <- df[df$annotation!='coding',]

for(i in 1:nrow(cluster_names)) {
  df$tissue[grep(cluster_names$cluster[i], df$annotation, ignore.case=T)] <- cluster_names$new_name[i]
}
df$annotation_type <- NA
df$annotation_type[grep('ATAC', df$annotation)] <- "ATAC peak"
df$annotation_type[grep('enhancer', df$annotation)] <- "Enhancer"
df$annotation_type[grep('TSS', df$annotation)] <- "Active TSS"

snp_counts <- annotations %>%
  dplyr::select(chrom, pos, locus, coding) %>%
  unique() %>%
  dplyr::group_by(locus) %>%
  dplyr::summarize(number_snps=n(),
                   coding_snps=sum(coding)) %>%
  dplyr::ungroup()


annotation_overlaps <- df %>%
  dplyr::group_by(locus, tissue, annotation_type) %>%
  dplyr::summarize(overlaps=sum(overlap)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(label=glue('{annotation_type} = {overlaps}')) %>%
  dplyr::group_by(locus, tissue) %>%
  dplyr::summarize(info=paste(sort(unique(label)), collapse='; ')) %>%
  dplyr::ungroup() %>%
  tidyr::spread(key=tissue, value=info)

tbl <- left_join(snp_counts, annotation_overlaps)

write.table(tbl, file='diamante-overlap-summary.csv', append = F, quote = T, sep = ',', row.names = F, col.names = T)
