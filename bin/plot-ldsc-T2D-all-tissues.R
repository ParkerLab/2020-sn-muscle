#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(glue)

args <- commandArgs(T)
CLUSTER_NAMES <- args[1]
LDSC_RESULT_FILES <- args[2:length(args)]

cluster_names <- read.table(CLUSTER_NAMES, head=T, sep='\t')
cluster_names$cluster <- paste0('cluster_', cluster_names$old_name)

load_result_file <- function(f) {
  tmp <- read.table(f, head = T, as.is = T, sep = '\t')
  NAME_RE <- '^tissues\\.(.*)\\.results$'
  tmp$trait <- gsub(NAME_RE, '\\1', basename(f))
  return(tmp)
}

results <- bind_rows(lapply(LDSC_RESULT_FILES, load_result_file))
results <- results[grep('L2_1', results$Category),]
results$Category <- gsub('L2_1', '', results$Category)
results$is_muscle_cell_type <- ifelse(grepl('cluster', results$Category), T, F)

for(i in 1:nrow(cluster_names)) {
  results$Category <- gsub(cluster_names$cluster[i], cluster_names$new_name[i], results$Category)
}


p <- ggplot(results) +
  geom_point(aes(x=Category, y=Coefficient, color=is_muscle_cell_type)) +
  geom_errorbar(aes(x=Category, ymin=Coefficient-1.96*Coefficient_std_error, ymax=Coefficient+1.96*Coefficient_std_error, color=is_muscle_cell_type), stat='identity', position='dodge') +
  theme_bw() +
  coord_flip() +
  facet_wrap(~trait, scales='free') +
  geom_hline(yintercept = 0, linetype='dashed') +
  ylab('Coefficient')
pdf('coefficients.pdf', width=25, height=25)
p
dev.off()
