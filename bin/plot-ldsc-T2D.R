#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(glue)

args <- commandArgs(T)
CLUSTER_NAMES <- args[1]
LDSC_RESULT_FILES <- args[2:length(args)]

# for testing
# LDSC_RESULT_FILES <- list.files('/lab/work/porchard/sn-muscle-project/work/ldsc-t2d-new-baseline-islets/results/partitioned-heritability', pattern='cluster.*.results', full.names = T)


cluster_names <- read.table(CLUSTER_NAMES, head=T, sep='\t')
cluster_names$cluster <- paste0('cluster_', cluster_names$old_name)

load_result_file <- function(f) {
  tmp <- read.table(f, head = T, as.is = T, sep = '\t')
  NAME_RE <- '^(cluster_\\d+)\\.(.*)\\.results$'
  tmp$trait <- gsub(NAME_RE, '\\2', basename(f))
  tmp$cluster <- gsub(NAME_RE, '\\1', basename(f))
  return(tmp)
}

results <- bind_rows(lapply(LDSC_RESULT_FILES, load_result_file))
results <- results[grep('L2_1', results$Category),]
results$Category <- gsub('L2_1', '', results$Category)
results <- left_join(results, cluster_names %>% dplyr::select(cluster, new_name) %>% dplyr::rename(cluster_name=new_name))

OTHER_CATEGORIES <- unique(grep('cluster_', results$Category, invert = T, value = T))


# plot the coefficient
p <- ggplot(results[!results$Category %in% OTHER_CATEGORIES,]) +
  geom_point(aes(x=cluster_name, y=Coefficient)) +
  geom_errorbar(aes(x=cluster_name, ymin=Coefficient-1.96*Coefficient_std_error, ymax=Coefficient+1.96*Coefficient_std_error), stat='identity', position='dodge') +
  theme_bw() +
  coord_flip() +
  facet_wrap(~trait, scales='free') +
  geom_hline(yintercept = 0, linetype='dashed') +
  ylab('Coefficient')
pdf('coefficient-muscle-cell-types.pdf', width=20, height=15)
p
dev.off()

for(CATEGORY in OTHER_CATEGORIES) {
  category_name <- gsub(' ', '_', CATEGORY)
  p <- ggplot(results[results$Category==CATEGORY,]) +
    geom_point(aes(x=cluster_name, y=Coefficient)) +
    geom_errorbar(aes(x=cluster_name, ymin=Coefficient-1.96*Coefficient_std_error, ymax=Coefficient+1.96*Coefficient_std_error), stat='identity', position='dodge') +
    theme_bw() +
    coord_flip() +
    facet_wrap(~trait, scales='free') +
    geom_hline(yintercept = 0, linetype='dashed') +
    ylab('Coefficient') +
    xlab('Other cell type in model')
  pdf(glue('coefficient-{category_name}.pdf'), width=20, height=15)
  print(p)
  dev.off()
}
