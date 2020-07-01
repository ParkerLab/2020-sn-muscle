#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(d3heatmap)
library(htmlwidgets)
library(glue)
library(ggplot2)

args <- commandArgs(T)
MODEL_RESULTS_HUMAN <- args[1]
MODEL_RESULTS_RAT <- args[2]
CLUSTER_NAMES <- args[3]
CELL_TYPE_DECODING <- args[4] # roadmap_cell_types.txt

DROP_TRAITS <- c()

models_human <- read.table(MODEL_RESULTS_HUMAN, head = T, sep = '\t')
models_human$species <- 'human'
models_rat <- read.table(MODEL_RESULTS_RAT, head = T, sep = '\t')
models_rat$species <- 'rat'
models <- bind_rows(models_human, models_rat)
models$cluster <- as.character(models$cluster)
models$cell_type <- as.character(models$cell_type)
models <- models[!models$cell_type %in% DROP_TRAITS,]

tissue_mappings <- read.table(CELL_TYPE_DECODING, head = T, as.is = T, sep = '\t', comment.char = '') %>%
  dplyr::rename(tissue_color=color, cell_type=eid, tissue=standardized_name) %>% 
  dplyr::select(tissue_color, cell_type, tissue)

# re-scale coefficients to between 0 and 1
for(cluster in unique(models$cluster)) {
  for(species in unique(models$species)) {
    models$coef[models$cluster==cluster & models$species==species] <- models$coef[models$cluster==cluster & models$species==species] / max(models$coef[models$cluster==cluster & models$species==species])
  }
}

#models$coef[models$coef<0] <- 0

models <- left_join(models, tissue_mappings %>% dplyr::select(cell_type, tissue))

# make the heatmap figure
cluster_names <- read.table(CLUSTER_NAMES, head=T, sep='\t') %>%
  mutate(cluster=glue('cluster_{old_name}'))
models <- left_join(models, cluster_names) %>%
  dplyr::select(new_name, cell_type, pvalue, coef, tissue, species) %>%
  rename(cluster=new_name)

models$cell_type <- as.character(models$cell_type)
models$tissue <- as.character(models$tissue)

KEEP <- c('Psoas muscle' = 'E100','MSC derived adipocytes' = 'E023', 'HUVEC cells (Endothelial)' = 'E122', 'Smooth muscle' = 'E111', 'Monocytes' = 'E029', 'Fetal trunk muscle' = 'E089')

for_heatmap_figure <- models[models$cell_type %in% KEEP,]
for_heatmap_figure$tissue <- sapply(for_heatmap_figure$cell_type, function(x){names(KEEP[KEEP==x])})
for_heatmap_figure$cluster <- factor(for_heatmap_figure$cluster, levels=rev(cluster_names$new_name), ordered=T)
for_heatmap_figure$tissue <- factor(for_heatmap_figure$tissue, levels=names(KEEP), ordered=T)
for_heatmap_figure <- for_heatmap_figure[order(for_heatmap_figure$cluster, for_heatmap_figure$tissue),]

p <- ggplot(for_heatmap_figure) +
  geom_tile(aes(x=cluster, y=tissue, fill=coef)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=55, hjust=1), legend.position = 'left') +
  scale_fill_gradient2(low='blue', mid='white', high='red', midpoint = mean(c(max(for_heatmap_figure$coef), min(for_heatmap_figure$coef)))) +
  guides(fill=guide_colorbar(title='Enhancer\nsimilarity\nscore')) +
  ylab('Roadmap cell type') +
  xlab('') +
  facet_wrap(~species, ncol=1) +
  coord_flip()
pdf('enhancer-similarity-heatmap.pdf', height=4, width=4)
p
dev.off()
