#!/usr/bin/env Rscript
library(glue)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)

# load in the UMAP
args <- commandArgs(T)
OUT <- args[1]
UMAP <- args[2]
TPM <- args[3]

# for testing
#UMAP <- '/lab/work/porchard/sn-muscle-project/work/liger/work/5e/69436debd8e2b48012b57bd0b851aa/factorization_k_20___factorization_lambda_50___norm_knnk_10___norm_resolution_0.4___umap_neighbors_10.dim.txt'
#TPM <- '/lab/work/porchard/sn-muscle-project/work/liger/work/5e/69436debd8e2b48012b57bd0b851aa/tpms.txt'

tpm <- read.table(TPM, head = F, as.is = T, sep = '\t', col.names = c('library', 'barcode', 'feature', 'score'))
umap <- read.table(UMAP, head = F, as.is = T, sep = '\t', col.names = c('library', 'barcode', 'dim1', 'dim2'))
umap <- umap[umap$library %in% tpm$library,]
umap <- umap %>%
  dplyr::mutate(nucleus=glue('{library}-{barcode}')) %>%
  dplyr::select(-library, -barcode)

tpm <- tpm %>%
	dplyr::mutate(nucleus=glue('{library}-{barcode}')) %>%
	dplyr::select(nucleus, feature, score) %>%
	tidyr::spread(key = feature, value=score, fill = 0) %>%
	tidyr::gather(key = feature, value=value, -nucleus)

NUMBER_MARKER_GENES <- length(unique(tpm$feature))

umap_with_genes <- left_join(umap, tpm)

if (!all(umap_with_genes$nucleus %in% tpm$nucleus)){
	warning('Dropping some nuclei without gene info')
}

umap_with_genes <- umap_with_genes[umap_with_genes$nucleus %in% tpm$nucleus,]

plots <- list()

for(i in unique(umap_with_genes$feature)) {
  p <- ggplot(umap_with_genes[umap_with_genes$feature==i,]) +
    theme_bw() +
    geom_point(aes(x = dim1, y = dim2, color = log(value+1)), alpha = 0.3, stroke = 0, size=1) +
    facet_wrap(~feature) +
    scale_color_gradient(low='grey', high='red', na.value = 'black') +
    xlab('UMAP Dim. 1') +
    ylab('UMAP Dim. 2') +
    guides(color=guide_colorbar(title='log(TPM+1)'))
  png(filename = gsub('.png', glue('-{i}.png'), OUT), height = 2, width = 3, res = 300, units = 'in')
  print(p)
  dev.off()
  plots[[length(plots)+1]] <- p
}

png(filename = OUT, height = 2*sqrt(NUMBER_MARKER_GENES), width = 4*sqrt(NUMBER_MARKER_GENES), res = 300, units = 'in')
cowplot::plot_grid(plotlist=plots)
dev.off()
