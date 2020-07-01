#!/usr/bin/Rscript

library(dplyr)
library(tidyr)

args <- commandArgs(T)
ld_file <- args[1]
pics_file <- args[2]
bdg_out <- args[3]

ld <- read.table(ld_file, head=F, as.is = T)[,1:3]
colnames(ld) <- c('chrom', 'pos', 'rsid')
ld$chrom <- paste0('chr', ld$chrom)
ld$start <- ld$pos-1
ld$end <- ld$pos

pics <- read.csv(pics_file)
pics <- pics[,c('Linked_SNP', 'PICS_probability')]
colnames(pics) <- c('rsid', 'ppa')

out <- left_join(ld, pics) %>%
  dplyr::select(chrom, start, end, ppa) %>%
  unique()
out <- out[order(out$start),]
write.table(out, file = bdg_out, append = F, quote = F, sep = '\t', row.names = F, col.names = F)
