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
all$library[grep('125589-', all$library)] <- '125589'

fragment_stats <- all %>%
  dplyr::group_by(library) %>%
  dplyr::summarize(median_fragments=round(median(fragments)),
                   mean_fragments=round(mean(fragments)))

sample_nuclei_counts <- all %>%
  dplyr::group_by(library, sample) %>%
  dplyr::summarize(number_nuclei=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sample=glue('nuclei_from_{sample}'))

library_info <- read.table(LIBRARY_LABELS, head=T, sep='\t', stringsAsFactors = F)
library_info$library <- gsub('125589-(hg19|rn6)', '125589', library_info$library)
library_info <- library_info %>%
  dplyr::select(library, fans_status, modality, experiment, loading_concentration) %>% unique() %>%
  left_join(sample_nuclei_counts) %>%
  left_join(fragment_stats)

DROP_LIBRARIES <- library_info$library[library_info$fans_status=='FANS' & library_info$modality=='ATAC']

table_out <- library_info %>% dplyr::filter(!library %in% DROP_LIBRARIES) %>%
  dplyr::select(-fans_status, -experiment) %>%
  tidyr::spread(key=sample, value=number_nuclei, fill=0)
table_out$library <- dplyr::recode(table_out$library,
				   '125589'='Human-rat mix (ATAC)',
				   '133151-hg19'='no FANS rep1 (ATAC)',
				   '133153-hg19'='no FANS rep2 (ATAC)',
				   '133155-hg19'='no FANS rep1 (RNA)',
				   '133156-hg19'='FANS rep1 (RNA)',
				   '133157-hg19'='no FANS rep2 (RNA)',
				   '133158-hg19'='FANS rep2 (RNA)',
				   '63_20_rna-hg19'='20k nuclei (RNA)',
				   '63_40_rna-hg19'='40k nuclei (RNA)',
				   '63_20-hg19'='20k nuclei (ATAC)',
				   '63_40-hg19'='40k nuclei (ATAC)')
colnames(table_out) <- dplyr::recode(colnames(table_out), 'nuclei_from_KSM1'='nuclei_from_HSM1', 'nuclei_from_KSM2'='nuclei_from_HSM2', 'nuclei_from_rat1'='nuclei_from_rat')
colnames(table_out) <- gsub('_', ' ', colnames(table_out))

write.table(table_out, file = 'nuclei-summary-stats.csv', append = F, quote = F, sep = ',', row.names = F, col.names = T)
