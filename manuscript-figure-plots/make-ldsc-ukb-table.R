#!/usr/bin/env Rscript
library(dplyr)
library(glue)

args <- commandArgs(T)
CLUSTER_NAMES <- args[1]
PHENOTYPE_TSV <- args[2]
LDSC_RESULT_FILES <- args[3:length(args)]

cluster_names <- read.table(CLUSTER_NAMES, head=T, sep='\t')
cluster_names$cluster <- paste0('cluster_', cluster_names$old_name)
cluster_names$new_name <- factor(cluster_names$new_name, levels=cluster_names$new_name, ordered=T)

load_result_file <- function(f) {
  tmp <- read.table(f, head = T, as.is = T, sep = '\t')
  NAME_RE <- '^.*\\.(.*)\\.results$'
  tmp$trait <- gsub(NAME_RE, '\\1', basename(f))
  return(tmp)
}

results <- bind_rows(lapply(LDSC_RESULT_FILES, load_result_file))

# for each trait, make one plot showing all annotations, baseline and custom
# for annotation, make one plot showing all traits

phenotypes <- read.table(gzfile(PHENOTYPE_TSV), head=T, as.is=T, sep='\t', quote='') %>%
  dplyr::select(phenotype, description) %>%
  unique() %>%
  rename(trait=phenotype)
results <- left_join(results, phenotypes)
results$Category <- gsub('L2_.*$', '', results$Category)
results <- results[grep('cluster', results$Category),]

results <- left_join(results, cluster_names %>% dplyr::select(cluster, new_name) %>% dplyr::rename(Category=cluster)) %>% dplyr::mutate(Category=new_name) %>% dplyr::select(-new_name)
results$description <- paste0(results$description, ' (', results$trait, ')')

# show z-score of all traits...
all <- results %>%
  dplyr::select(Category, Coefficient_z.score, description) %>%
  tidyr::spread(key=Category, value=Coefficient_z.score)

write.table(all, file = 'LDSC-UKB-Z-scores.tsv', append = F, quote = F, sep = '\t', row.names = F, col.names = T)
