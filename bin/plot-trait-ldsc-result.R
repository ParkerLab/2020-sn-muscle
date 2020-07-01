#!/usr/bin/env Rscript
library(tidyr)
library(dplyr)
library(ggplot2)
library(glue)

args <- commandArgs(T)
LDSC_RESULT_FILE <- args[1]
CLUSTER_NAMES <- args[2]

cluster_names <- read.table(CLUSTER_NAMES, head=T, as.is=T, sep='\t')

load_result_file <- function(f) {
  tmp <- read.table(f, head = T, as.is = T, sep = '\t')
  NAME_RE <- '^tissues\\.(.*)\\.results$'
  tmp$trait <- gsub(NAME_RE, '\\1', basename(f))
  return(tmp)
}

results <- load_result_file(LDSC_RESULT_FILE)
results$col <- 'LDSC baseline'
results$col[grep('L2_1', results$Category)] <- 'Other tissue open chromatin'
results$col[grep('cluster', results$Category)] <- 'Muscle cell type open chromatin'
results$Category <- gsub('L2_[01]', '', results$Category)

for(i in 1:nrow(cluster_names)) {
  cl <- cluster_names$old_name[i]
  results$Category[results$Category==glue('cluster_{cl}')] <- cluster_names$new_name[i]
}

# read in the phenotype codes
phenotypes <- read.table('/lab/work/porchard/sn-muscle-project/data/ukb-summary-stats/manifest.tsv', head=T, as.is = T, sep = '\t', comment.char='', quote='')
phenotypes <- phenotypes[!is.na(phenotypes$Phenotype.Code),]
phenotypes <- phenotypes[phenotypes$Sex=='both_sexes',]
phenotypes$trait <- phenotypes$Phenotype.Code
results <- left_join(results, phenotypes)
results <- dplyr::select(results, -wget.command, -Dropbox.File, -Phenotype.Code, -UK.Biobank.Data.Showcase.Link, -File)
trait_name <- gsub(' ', '_', unique(results$Phenotype.Description))
trait_name <- gsub('/', '_', trait_name)
trait_name <- gsub('%', '_', trait_name)
trait_name <- gsub(',', '', trait_name)
trait_name <- gsub('\\(', '', trait_name)
trait_name <- gsub('\\)', '', trait_name)

results <- results[order(results$Coefficient_z.score),]
results$Category <- factor(results$Category, levels=results$Category, ordered=T)

p <- ggplot(results) +
  geom_errorbar(aes(x = Category, color=col, ymin=Coefficient-1.96*Coefficient_std_error, ymax=Coefficient+1.96*Coefficient_std_error)) +
  geom_point(aes(x = Category, y=Coefficient, color=col)) +
  theme_bw() +
  xlab('Annotation') + ylab('Coefficient') +
  coord_flip() +
  geom_hline(yintercept = 0, linetype='dashed')
pdf(glue('{trait_name}.pdf'), width=8, height=15)
p
dev.off()



