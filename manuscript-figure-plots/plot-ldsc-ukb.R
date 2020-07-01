#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(glue)
library(htmlwidgets)
library(ComplexHeatmap)
# library(d3heatmap)
library(gplots)

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

phenotypes <- read.table(gzfile(PHENOTYPE_TSV), head=T, as.is=T, sep='\t', quote='') %>%
  dplyr::select(phenotype, description) %>%
  unique() %>%
  dplyr::rename(trait=phenotype)
results <- left_join(results, phenotypes)

# for each tissue, get the traits where it ranks #1
results <- results[grep('L2_1', results$Category),]
results$Category <- gsub('L2_.*$', '', results$Category)
top_tissue_per_trait <- results %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(max_z_score=max(Coefficient_z.score)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(Coefficient_z.score==max_z_score)

TOP_IN_ONE_OF_OUR_CELL_TYPES <- top_tissue_per_trait$trait[grepl('cluster', top_tissue_per_trait$Category) & top_tissue_per_trait$Coefficient_z.score>=3]

results <- results[grep('cluster', results$Category),]
results <- left_join(results, cluster_names %>% dplyr::select(cluster, new_name) %>% dplyr::rename(Category=cluster)) %>% dplyr::mutate(Category=new_name) %>% dplyr::select(-new_name)
results$pval <- sapply(pnorm(results$Coefficient_z.score, lower.tail = F)*2, function(x){min(x, 1)})
NUMBER_TRAITS <- length(unique(results$trait))
BONFERRONI_THRESHOLD <- 0.05/nrow(results)
results$significant <- ifelse(results$pval <= BONFERRONI_THRESHOLD, 'Sign. after bonferroni', 'N.S.')

# plot all things that are top in one of our cell types
top_in_our_cells <- results[results$trait %in% TOP_IN_ONE_OF_OUR_CELL_TYPES,]
top_in_our_cells <- top_in_our_cells %>%
  dplyr::select(Category, Coefficient_z.score, description) %>%
  tidyr::spread(key=Category, value=Coefficient_z.score)
rownames(top_in_our_cells) <- top_in_our_cells$description
top_in_our_cells <- top_in_our_cells %>% dplyr::select(-description) %>% as.matrix()
rownames(top_in_our_cells) <- sapply(rownames(top_in_our_cells), function(x){paste(strwrap(x, width=75), collapse='\n')})
pdf('UKB-top-in-our-cell-types.pdf', height=13, width = 7)
hm <- Heatmap(top_in_our_cells, name = 'LDSC Z-score', heatmap_legend_param = list(direction = "horizontal"), show_row_dend = F, column_names_max_height = unit(10, 'in'), show_column_dend = F, width=unit(1.5, 'in'), heatmap_height = unit(11, "in")) #, heatmap_height = unit(5, "in"))
draw(hm, heatmap_legend_side = "top")
dev.off()

# plot all things that are sign. after BY correction
results$significant <- ifelse(p.adjust(results$pval, method='BY') <= 0.05, 'Sign. after BY', 'N.S.')
by_traits <- unique(results$trait[results$significant=='Sign. after BY'])
by_sign <- results[results$trait %in% by_traits,]
by_sign$description[by_sign$description=='Blood clot, DVT, bronchitis, emphysema, asthma, rhinitis, eczema, allergy diagnosed by doctor: Hayfever, allergic rhinitis or eczema'] <- 'Diagnosed by doctor: Hayfever, allergic rhinitis or eczema'

by_sign <- by_sign %>%
  dplyr::select(Category, Coefficient_z.score, description) %>%
  tidyr::spread(key=Category, value=Coefficient_z.score)
rownames(by_sign) <- by_sign$description
by_sign <- by_sign %>% dplyr::select(-description) %>% as.matrix()
pdf('UKB-LDSC-by.pdf', height=7, width = 11)
hm <- Heatmap(t(by_sign), cluster_rows = F, name = 'LDSC Z-score', heatmap_legend_param = list(direction = "horizontal"), height=unit(1.3, 'in'), column_names_max_height = unit(10, 'in'), show_row_dend = F, show_column_dend = F, heatmap_width = unit(10, "in"))
draw(hm, heatmap_legend_side = "top")
dev.off()
