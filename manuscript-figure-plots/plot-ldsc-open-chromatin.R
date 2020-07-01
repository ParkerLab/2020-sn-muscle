#!/usr/bin/env Rscript
library(tidyr)
library(dplyr)
library(ggplot2)
library(glue)

args <- commandArgs(T)
LDSC_RESULT_FILE <- args[1]
CLUSTER_NAMES <- args[2]
OUT <- args[3]
TITLE <- args[4]

capitalize <- function(s) {
  ### capitalize the first letter
  tmp <- strsplit(s, '')[[1]]
  tmp[1] <- toupper(tmp[1])
  tmp <- paste(tmp, collapse='')
  return(tmp)
}

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

results <- results[order(results$Coefficient_z.score),]
results <- results[results$col!='LDSC baseline',]
results$Category <- gsub('_', ' ', results$Category)
results$Category <- tolower(results$Category)
results$Category <- sapply(results$Category, capitalize)
results$Category[results$Category=='Beta atac'] <- 'Beta cells'
results$Category <- gsub(' i ', ' I ', results$Category)
results$Category <- gsub(' ii ', ' II ', results$Category)
results$Category <- factor(results$Category, levels=results$Category, ordered=T)
results <- results[order(results$Category),]
COLORS <- ifelse(results$col=='Muscle cell type open chromatin', 'red', 'black')

p <- ggplot(results) +
  #geom_errorbar(aes(x = Category, ymin=Coefficient-1.96*Coefficient_std_error, ymax=Coefficient+1.96*Coefficient_std_error)) +
  geom_linerange(aes(x = Category, ymin=Coefficient-1.96*Coefficient_std_error, ymax=Coefficient+1.96*Coefficient_std_error)) +
  geom_point(aes(x = Category, y=Coefficient)) +
  theme_bw() +
  xlab('') + ylab('LDSC coefficient') +
  coord_flip() +
  geom_hline(yintercept = 0, linetype='dashed') +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust=0.5), axis.text.y = element_text(colour = COLORS)) +
  ggtitle(TITLE)
pdf(OUT, width=4, height=6)
p
dev.off()



