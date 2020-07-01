#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(glue)

args <- commandArgs(T)
LIBRARY_LABELS <- args[1]
NUCLEI_TO_INDIVIDUAL <- args[2]
RNA_METRICS <- args[3]
ATAC_METRICS <- args[4]


rna_fragment_counts <- read.table(RNA_METRICS, sep = '\t', head=F, stringsAsFactors = F, col.names = c('library', 'barcode', 'fragments'))
atac_fragment_counts <- read.table(ATAC_METRICS, sep='\t', head=F, stringsAsFactors = F, col.names = c('nucleus', 'metric', 'value')) %>%
  dplyr::filter(metric=='hqaa') %>%
  dplyr::mutate(value=as.numeric(value),
                library=gsub('(.*)-(.*)', '\\1', nucleus),
                barcode=gsub('(.*)-(.*)', '\\2', nucleus),
                fragments=value/2) %>%
  dplyr::select(-metric, -nucleus, -value)

individuals <- read.table(NUCLEI_TO_INDIVIDUAL, head=F, stringsAsFactors = F, sep='\t', col.names = c('library', 'barcode', 'sample'))
all <- left_join(individuals, bind_rows(atac_fragment_counts, rna_fragment_counts))

library_info <- read.table(LIBRARY_LABELS, head=T, sep='\t', stringsAsFactors = F)
library_info <- library_info %>%
  dplyr::select(library, fans_status, modality, experiment, loading_concentration) %>% unique()

all <- left_join(all, library_info)
all$genome <- gsub('.*-', '', all$library)
all <- all[!(all$modality=='ATAC' & all$fans_status=='FANS'),]
per_species_per_modality <- all[,c('library', 'barcode', 'genome', 'modality', 'fragments')] %>%
  unique() %>%
  dplyr::group_by(genome, modality) %>%
  dplyr::summarize(median_fragments=round(median(fragments)),
                   mean_fragments=round(mean(fragments))) %>%
  dplyr::ungroup()

write.table(per_species_per_modality, file = 'per-species-and-modality-fragment-stats.csv', append = F, quote = F, sep = ',', row.names = F, col.names = T)
