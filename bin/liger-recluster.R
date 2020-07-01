#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option(c("--rda"), action = "store", type = "character", default = "", help = "[Required] Path to input .Rda file"),
  make_option(c("--cluster_file"), action = "store", type = "character",help = "[Required] Path to output cluster file"),
  make_option(c("--resolution"), action = "store", default=0.1, type = "numeric",help = "[Required] Resolution to use for the re-clustering (default: 0.1)"),
  make_option(c("--k"), action = "store", default=0.1, type = "numeric",help = "[Required] K nearest neighbors for the re-clustering (default: 20)"),
  make_option(c("--out"), action = "store", type = "character", default = '', help = "[Required] Output .Rda file.")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts <- parse_args(option_parser)

library(liger)
library(glue)
library(Seurat)

load(opts$rda)

# re-cluster using clusterLouvainJaccard
# need seurat v2 here
int.muscle <- clusterLouvainJaccard(int.muscle, k.param=opts$k, resolution=opts$resolution)

df <- data.frame(cluster=int.muscle@clusters)
df$library <- gsub('(.*)-.*', '\\1', rownames(df))
df$barcode <- gsub('.*-', '', rownames(df))

write.table(df[,c('library', 'barcode', 'cluster')], file=opts$cluster_file, append=F, quote=F, sep='\t', row.names=F, col.names=F)

save(int.muscle, file=opts$out)
