#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(glue)

args <- commandArgs(T)
CLUSTER_NAMES <- args[1]
LDSC_RESULT_FILES <- args[2:length(args)]

# for testing
# LDSC_RESULT_FILES <- list.files('/lab/work/porchard/sn-muscle-project/work/ldsc-ukb-new-baseline/results/partitioned-heritability', pattern='.results', full.names = T)

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

# split phenotypes by category???


# split the traits by category too...?
fields <- read.table('/lab/work/porchard/sn-muscle-project/data/ukb-schemas/field.txt', head=T, as.is=T, sep='\t', quote='', comment.char = '')
fields$field_id <- as.character(fields$field_id)
categories <- read.table('/lab/work/porchard/sn-muscle-project/data/ukb-schemas/category.txt', head=T, as.is=T, sep='\t', quote='', comment.char = '')

phenotypes <- read.table(gzfile('/lab/work/porchard/daily/2019-12-28/confident-phenotypes/ukb31063_h2_all.02Oct2019.tsv.gz'), head=T, as.is=T, sep='\t', quote='') %>%
  dplyr::select(phenotype, description) %>%
  unique() %>%
  rename(trait=phenotype)
phenotypes$field <- gsub('_.*', '', phenotypes$trait)
stopifnot(all(phenotypes$field %in% fields$field_id))
phenotypes <- left_join(phenotypes, fields %>% dplyr::select(field_id, title, main_category) %>% rename(field=field_id))


results <- left_join(results, phenotypes)
results$Category <- gsub('L2_.*$', '', results$Category)
results <- results[results$Category!='base',]

for(i in 1:nrow(cluster_names)) {
  results$Category <- gsub(cluster_names$cluster[i], cluster_names$new_name[i], results$Category)
}

# for each annotation, make one plot showing all traits
for(ANNOTATION in unique(results$Category)) {
  annotation_name <- gsub(' ', '_', ANNOTATION)
  tmp <- results[results$Category==ANNOTATION,]
  tmp <- tmp[rev(order(tmp$Coefficient_z.score)),]
  tmp$description <- paste0(tmp$description, ' (', tmp$trait, ')')
  tmp$description <- factor(tmp$description, levels=tmp$description, ordered=T)
  p <- ggplot(tmp) +
    geom_point(aes(x=description, y=Coefficient)) +
    geom_errorbar(aes(x=description, ymin=Coefficient-1.96*Coefficient_std_error, ymax=Coefficient+1.96*Coefficient_std_error), stat='identity', position='dodge') +
    theme_bw() +
    coord_flip() +
    geom_hline(yintercept = 0, linetype='dashed') +
    ylab('Coefficient')
  pdf(glue('all-traits-{annotation_name}.pdf'), width=12, height=35)
  print(p)
  dev.off()
  
  p <- ggplot(tmp[tmp$Coefficient/tmp$Coefficient_std_error>=1.96,]) +
    geom_point(aes(x=description, y=Coefficient)) +
    geom_errorbar(aes(x=description, ymin=Coefficient-1.96*Coefficient_std_error, ymax=Coefficient+1.96*Coefficient_std_error), stat='identity', position='dodge') +
    theme_bw() +
    coord_flip() +
    geom_hline(yintercept = 0, linetype='dashed') +
    ylab('Coefficient')
  pdf(glue('all-significant-traits-{annotation_name}.pdf'), width=12, height=10)
  print(p)
  dev.off()
}

# for each trait, make one plot showing all annotations, baseline and custom
for(TRAIT in unique(results$trait)) {
  tmp <- results[results$trait==TRAIT,]
  trait_name <- unique(paste(tmp$description, tmp$trait))
  trait_name <- gsub("'", '', trait_name)
  trait_name <- gsub(" ", '_', trait_name)
  trait_name <- gsub("\\(", '', trait_name)
  trait_name <- gsub("\\)", '', trait_name)
  trait_name <- gsub(":", '', trait_name)
  trait_name <- gsub("\\/", '', trait_name)
  tmp <- tmp[order(tmp$Coefficient_z.score),]
  tmp$Category <- factor(tmp$Category, levels=tmp$Category, ordered=T)
  p <- ggplot(tmp) +
    geom_point(aes(x=Category, y=Coefficient)) +
    geom_errorbar(aes(x=Category, ymin=Coefficient-1.96*Coefficient_std_error, ymax=Coefficient+1.96*Coefficient_std_error), stat='identity', position='dodge') +
    theme_bw() +
    coord_flip() +
    facet_wrap(~description, scales='free') +
    geom_hline(yintercept = 0, linetype='dashed') +
    ylab('Coefficient')
  pdf(glue('all-annotations-{trait_name}.pdf'), width=8, height=10)
  print(p)
  dev.off()
}


# for our custom annotations (and/or all DHS annotations as well?) plot out significant traits...
# could sort traits based on what they are significant in....?


