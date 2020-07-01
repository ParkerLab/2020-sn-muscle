#!/usr/bin/env Rscript

args <- commandArgs(T)
CONFIDENCES <- args[1]
MANIFEST <- args[2]

confidences <- read.table(gzfile(CONFIDENCES), head=T, as.is=T, sep='\t', quote='')
confidences <- confidences[,grep('Prop_h2', colnames(confidences), invert = T)]
confidences <- confidences[,grep('\\.Enrichment', colnames(confidences), invert = T)]
confidences <- confidences[,grep('\\.Coefficient', colnames(confidences), invert = T)]
confidences <- confidences[,grep('Prop_SNPs', colnames(confidences), invert = T)]
KEEP <- confidences$phenotype[confidences$confidence %in% c('high') & confidences$h2_z >= 7 & confidences$h2_observed>=0.01 & confidences$keep==T]

manifest <- read.table(MANIFEST, head=T, as.is=T, sep='\t', quote='', comment.char = '')
manifest <- manifest[manifest$Phenotype.Code!='N/A',]
manifest <- manifest[manifest$Sex=='both_sexes',]
AVAILABLE <- manifest$Phenotype.Code

MISSING <- KEEP[!KEEP %in% AVAILABLE]

manifest_subset <- manifest[manifest$Phenotype.Code %in% KEEP,]
DUPLICATES <- table(manifest_subset$Phenotype.Code)
DUPLICATES <- DUPLICATES[DUPLICATES>1]
# manifest_subset[manifest_subset$Phenotype.Code %in% names(DUPLICATES),]

manifest_subset <- manifest_subset[manifest_subset$md5s!='<pending>',c('Phenotype.Code', 'Phenotype.Description', 'File', 'wget.command', 'md5s')]
write.table(manifest_subset, 'manifest-subset.tsv', col.names = F, append = F, quote = F, sep = '\t', row.names = F)

keep <- data.frame(p=KEEP[KEEP %in% AVAILABLE])
write.table(keep, 'ukb-traits-new.txt', col.names = F, append = F, quote = F, sep = '\t', row.names = F)
