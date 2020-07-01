#!/usr/bin/env Rscript
## This function will perform DoubletFinder analyses on the input dataset and
## output a text file with the barcodes flagged as doublets. In addition, it
## will generate Seurat plots for marker genes before and after removing
## doublets to facilitate QC.

options(stringsAsFactors = F)
seed <- 87532163

library(optparse)
library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)
library(Seurat) # must be v3...
library(DoubletFinder)
library(glue)


#source("~/src/R_rainclouds.R")
theme_set(theme_bw(base_size = 12))

option_list <- list(
  make_option(c("-i", "--input"), type = "character", 
              help = "Input feature file with RNA/ATAC counts"),
  make_option(c("-o", "--outdir"), type = "character",
              help = "Output directory"),
  make_option("--sample", type = "character", 
              help = "Sample name"),
  make_option("--max_dims", type = "character", default = 20, 
              help = paste0(
                "Number of data dimensions to account for in clustering ",
                "(higher numbers yield more separate clusters)"
              ))
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", 
                              option_list = option_list, add_help_option = T)
opts <- parse_args(option_parser)

# Read parameters
infile <- opts$input
outdir <- opts$outdir
sample <- opts$sample
max_dims <- opts$max_dims
## Genes used as markers
gene_list <- c("MYH1", "MYH7", "TNNT1", "TNNT3", "VWF", 
               "PDGFRA", "PAX7", "CD163", 'FBN1')

# # Test
#infile <- "/lab/work/porchard/sn-muscle-project/work/doubletfinder/results/counts/63_20-hg19.filtered-counts.txt"
#infile <- "counts.txt"
#sample <- "HPAP045"
#outdir <- "."
#max_dims <- 20

print(sprintf("Processing data with %s dimensions", max_dims))


counts <- read.table(infile, head = F, sep = '\t', col.names = c('library', 'barcode', 'gene', 'count'), colClasses = c('character', 'character', 'character', 'numeric'))
counts <- counts %>%
  dplyr::group_by(library, barcode, gene) %>%
  dplyr::summarize(count=sum(count)) %>%
  dplyr::ungroup()
counts_matrix <- counts %>%
	  dplyr::mutate(nucleus=glue('{library}-{barcode}')) %>%
	    dplyr::select(-library, -barcode) %>%
	      tidyr::spread(key = nucleus, value = count, fill = 0)

counts_matrix <- as.data.frame(counts_matrix[!is.na(counts_matrix$gene),])
rownames(counts_matrix) <- counts_matrix$gene
counts_matrix <- dplyr::select(counts_matrix, -gene)
counts_matrix <- as.matrix(counts_matrix)

counts_data <- counts_matrix
print('Making seurat object')
# Create and process Seurat object from sparse counts matrix
sobj <- CreateSeuratObject(counts_data, project = "test", 
                           min.cells = 3, min.features = 200)
sobj <- NormalizeData(sobj)
sobj <- ScaleData(sobj)
#sobj <- FindVariableGenes(sobj, top.genes = 2000)
sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)
#VARFEATURES <- VariableFeatures(sobj)
#VariableFeatures(sobj, assay = NULL, selection.method = NULL) <- VARFEATURES
sobj <- RunPCA(sobj)
# ElbowPlot(sobj)
sobj <- RunUMAP(sobj, dims = 1:max_dims)
sobj <- FindNeighbors(sobj, dims = 1:max_dims)
sobj <- FindClusters(sobj, resolution = 0.8)

# Run DoubletFinder
## Estimate DoubletFinder parameters from data
sweep.res.list <- paramSweep_v3(sobj, PCs = 1:max_dims, sct = FALSE, num.cores = 1)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
## Make a list with all the parameters
DF_pars <- list()
DF_pars$pN = 0.25
DF_pars$pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>% 
  mutate(pK = as.numeric(as.character(pK))) %>% 
  pull(pK)
DF_pars$nExp_poi <- round(0.075 * length(sobj@meta.data$orig.ident))
homotypic.prop <- modelHomotypic(sobj@meta.data$seurat_clusters)
DF_pars$nExp_poi.adj <- round(DF_pars$nExp_poi * (1 - homotypic.prop))
## DoubletFinder appends columns to the metadata named after the parameters
pANN_lab = paste("pANN", DF_pars$pN, DF_pars$pK, DF_pars$nExp_poi, sep = "_")
DF_col <- gsub("pANN", "DF.classifications", pANN_lab)
DF_col2 <- paste("DF.classifications", DF_pars$pN, DF_pars$pK, 
                 DF_pars$nExp_poi.adj, sep = "_")
## Flag doublets with two levels of stringency (high and lower, resp.)
sobj <- doubletFinder_v3(sobj, PCs = 1:max_dims, pN = DF_pars$pN, 
                         pK = DF_pars$pK, nExp = DF_pars$nExp_poi, 
                         reuse.pANN = FALSE, sct = FALSE)
sobj <- doubletFinder_v3(sobj, PCs = 1:max_dims, pN = DF_pars$pN, 
                         pK = DF_pars$pK, nExp = DF_pars$nExp_poi.adj, 
                         reuse.pANN = pANN_lab, sct = FALSE)


# Format data, collect results, and make plots
sobj@meta.data[,"df"] <- sobj@meta.data[,DF_col]
sobj@meta.data[,"df2"] <- sobj@meta.data[,DF_col2]

md <- sobj@meta.data
singlets <- rownames(md)[md$df == "Singlet"]
assignments <- data.frame(nucleus=rownames(md), assignment=md$df)
assignments$library <- gsub('(.*)-(.*)', '\\1', assignments$nucleus)
assignments$barcode <- gsub('(.*)-(.*)', '\\2', assignments$nucleus)
assignments$assignment <- as.character(assignments$assignment)
assignments$assignment[assignments$assignment=='Singlet'] <- 'singlet'
assignments$assignment[assignments$assignment=='Doublet'] <- 'doublet'
write.table(assignments[,c('library', 'barcode', 'assignment')], file="assignments.txt", quote = F, append = F, sep = '\t', row.names = F, col.names = F)

## Count flagged doublets
doublets_summary <- data.frame(
  sample = sample,
  homotypic.prop,
  pK = DF_pars$pK,
  high_stringency = sum(sobj@meta.data$df == "Doublet"),
  low_stringency = sum(sobj@meta.data$df2 == "Doublet")
)
outfile <- file.path(outdir, "doublet_summary.txt")
write.table(doublets_summary, outfile, sep = "\t", quote = F, row.names = F)

## Output barcodes not flagged as singlets
outfile <- file.path(outdir, "barcodes_singlets.txt")
singlets %>% gsub(".*_", "", .) %>% 
  write.table(outfile, col.names = F, row.names = F, quote = F)

## Plot UMAPs
outplot <- file.path(outdir, "qc_umaps.pdf")
cairo_pdf(outplot, 6, 4)
DimPlot(sobj, group.by="df", reduction="umap", pt.size=0.5) +
  facet_wrap(~ df) +
  ggtitle("nExp_poi")
FeaturePlot(sobj, features = "nCount_RNA") +
  facet_wrap(~ sobj@meta.data$df)
FeaturePlot(sobj, features = "nFeature_RNA") +
  facet_wrap(~ sobj@meta.data$df)
dev.off()

## Plot boxplots of gene/feature counts
# outplot <- file.path(outdir, "qc_counts.pdf")
# cairo_pdf(outplot, 3, 3)
# md %>% 
#   mutate(df = factor(df)) %>% 
#   ggplot(aes(x = df, y = nFeature_RNA, group = df)) + 
#   geom_flat_violin(position = position_nudge(x = .25, y = 0),
#                    adjust = 0.8, trim = FALSE, alpha = .5) +
#   geom_point(position = position_jitter(width = .15), size = .25, alpha = .05) +
#   geom_boxplot(aes(x = as.numeric(df) + 0.25),
#                outlier.shape = NA, alpha = 1, width = .2) +
#   theme()
# md %>% 
#   mutate(df = factor(df)) %>% 
#   ggplot(aes(x = df, y = nCount_RNA, group = df)) + 
#   geom_flat_violin(position = position_nudge(x = .25, y = 0),
#                    adjust = 0.8, trim = FALSE, alpha = .5) +
#   geom_point(position = position_jitter(width = .15), size = .25, alpha = .05) +
#   geom_boxplot(aes(x = as.numeric(df) + 0.25),
#                outlier.shape = NA, alpha = 1, width = .2) +
#   scale_y_continuous(trans = "log10", labels = comma) +
#   theme()
# dev.off()


# Recluster data without doublets
## Remove doublets and save output
counts_data_filt <- counts_data[,which(colnames(counts_data) %in% singlets)]
outfile <- file.path(outdir, basename(infile))  # {rna,atac}_processed.rds
saveRDS(counts_data_filt, file = outfile)

## Create and process Seurat object from filtered counts matrix
sobj_filt <- CreateSeuratObject(counts_data_filt, project = "test", 
                                min.cells = 3, min.features = 200)
sobj_filt <- NormalizeData(sobj_filt)
sobj_filt <- ScaleData(sobj_filt)
sobj_filt <- FindVariableFeatures(sobj_filt, selection.method = "vst", nfeatures = 2000)
sobj_filt <- RunPCA(sobj_filt)
sobj_filt <- RunUMAP(sobj_filt, dims = 1:max_dims)
sobj_filt <- FindNeighbors(sobj_filt, dims = 1:max_dims)
sobj_filt <- FindClusters(sobj_filt, resolution = 0.8)

## Plot clustering UMAPs
outplot <- file.path(outdir, "clustering.pdf")
cairo_pdf(outplot, 3, 3)
DimPlot(sobj, reduction = "umap") +
  ggtitle("Doublets retained")
DimPlot(sobj_filt, reduction = "umap") +
  ggtitle("Doublets removed")
dev.off()

## Plot marker genes on UMAP
outplot <- file.path(outdir, "markers_umap.pdf")
cairo_pdf(outplot, 9, 8)
FeaturePlot(sobj, features = gene_list)
FeaturePlot(sobj_filt, features = gene_list)
dev.off()

## Plot marker genes expression across clusters
outplot <- file.path(outdir, "markers_vln.pdf")
cairo_pdf(outplot, 9, 8)
VlnPlot(sobj, features = gene_list, slot = "counts", log = TRUE, pt.size = .1)
VlnPlot(sobj_filt, features = gene_list, slot = "counts", log = TRUE, pt.size = .1)
dev.off()
