#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option(c("--rda"), action = "store", type = "character", default = "", help = "[Required] Path to input .Rda file"),
  make_option(c("--norm_resolution"), action = "store", type = "numeric", default = 0.4, help = "[Required] resolution parameter for quantile normalization"),
  make_option(c("--norm_knnk"), action = "store", type = "numeric", default = 20, help = "[Required] knn_k parameter for quantile normalization"),
  make_option(c("--dim"), action = "store", type = "character", default = '', help = "[Required] Output dim file"),
  make_option(c("--out"), action = "store", type = "character", default = '', help = "[Required] Output .Rda file.")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts <- parse_args(option_parser)

# for testing
#opts <- list(
#   'mats' = paste(list.files('/lab/work/porchard/sn-muscle-project/work/downstream-stricter-qc/results/liger/round-1/input-merged/', full.names=T), collapse=','),
#   'factorization_k' = 20,
#   'factorization_lambda' = 5,
#   'norm_resolution'=0.4,
#   'norm_knnk' = 20,
#   'cluster_file'='test-cluster-file.txt',
#   'dim'='test-dim.txt',
#   'out'='test.Rda'
#)


library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(liger)
library(glue)

print('Reference dataset:')
SELECT_GENES <- 'RNA-KSM1.liger'
print(SELECT_GENES)

load(opts$rda)

# remove ribosomal factor(s) or other technical factors
markers <- getFactorMarkers(int.muscle, num.genes = 10)
shared <- markers$shared
shared$is_ribosomal <- F
shared$is_ribosomal[shared$is_ribosomal==F] <- grepl('^RPS', shared$gene[shared$is_ribosomal==F])
shared$is_ribosomal[shared$is_ribosomal==F] <- grepl('^RPL', shared$gene[shared$is_ribosomal==F])
RIBOSOMAL_FACTORS <- shared %>% dplyr::group_by(factor_num) %>% dplyr::summarize(number_ribosomal=sum(is_ribosomal)) %>% dplyr::filter(number_ribosomal>0)
KEEP_FACTORS <- 1:max(shared$factor_num)
DROP_FACTORS <- c(3, 5)
KEEP_FACTORS <- KEEP_FACTORS[!KEEP_FACTORS %in% DROP_FACTORS]
print('Using factors:')
print(KEEP_FACTORS)
shared$keep_factor <- shared$factor_num %in% KEEP_FACTORS
write.table(shared, file = 'shared-factor-markers.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = T)

# it's possible to have cells that only load on the factors that have been removed -- that will lead to errors in the downstream processing
# therefore, remove them
nonzero_cells <- c()
for(i in 1:length(int.muscle@H)) {
  nonzero_cells = c(nonzero_cells, names(which(rowSums(int.muscle@H[[i]][,KEEP_FACTORS]>0)>0)))
}

int.muscle = subsetLiger(int.muscle, cells.use = nonzero_cells)
int.muscle <- quantileAlignSNF(int.muscle, knn_k = opts$norm_knnk, center=T, resolution = opts$norm_resolution, small.clust.thresh = opts$norm_knnk, dims.use = KEEP_FACTORS, ref_dataset=SELECT_GENES) 
int.muscle <- runUMAP(int.muscle, dims.use = KEEP_FACTORS, n_neighbors=15)
umap <- as.data.frame(int.muscle@tsne.coords)
colnames(umap) <- c('dim_1', 'dim_2')
umap$nucleus <- rownames(umap)
umap$library <- gsub('(.*)-(.*)', '\\1', umap$nucleus)
umap$barcode <- gsub('(.*)-(.*)', '\\2', umap$nucleus)
umap <- dplyr::select(umap, library, barcode, dim_1, dim_2)
write.table(umap, file = opts$dim, append = F, quote = F, sep = '\t', row.names = F, col.names = F)

save(int.muscle, file=opts$out)
