#!/usr/bin/env Rscript
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggseqlogo)
library(glue)


args <- commandArgs(T)
EXPLAINED_FILES <- args


read_explain_file <- function(f) {
  FILENAME_RE <- '(cluster_.*)\\.(ref|alt).explained.reformatted.txt'
  tmp <- read.table(f, head=F, sep='\t')
  colnames(tmp) <- c('snp', 'pos', 'nuc', 'score')
  tmp$cluster <- gsub(FILENAME_RE, '\\1', basename(f))
  tmp$ref_alt <- gsub(FILENAME_RE, '\\2', basename(f))
  return(tmp)
}


explained <- bind_rows(lapply(EXPLAINED_FILES, read_explain_file))

# now plot, for each cluster and SNP, ref/alt importance scores above one another
#SNP <- 'chr5:53271420:G:A'
#CLUSTER <- 'cluster_2'

to_matrix <- function(cluster, snp, ref_or_alt, positions='all') {
  tmp <- explained[explained$cluster==cluster & explained$snp==snp & explained$ref_alt==ref_or_alt,]
  if (positions != 'all') {
    # keep the middle 'positions' positions
    middle <- ceiling(max(tmp$pos) / 2)
    lower <- middle - (0.5 * positions)
    upper <- middle + (0.5 * positions)
    tmp <- tmp[tmp$pos >= lower & tmp$pos <= upper,]
  }
  tmp <- tmp[order(tmp$pos, tmp$nuc),] %>%
    dplyr::select(pos, nuc, score) %>%
    tidyr::spread(key=pos, value=score) %>%
    dplyr::select(-nuc) %>%
    as.matrix()
  colnames(tmp) <- 1:ncol(tmp)
  rownames(tmp) <- c('A', 'C', 'G', 'T')
  return(tmp)
}

for(CLUSTER in unique(explained$cluster)) {
  for(SNP in unique(explained$snp)) {
    ref <- to_matrix(CLUSTER, SNP, 'ref', positions=50)
    alt <- to_matrix(CLUSTER, SNP, 'alt', positions=50)
    ref_minus_alt <- ref - alt
    p <- ggseqlogo(list('With ref. allele'=ref, 'With alt. allele'=alt, 'Ref - alt'=ref_minus_alt), method='custom', seq_type='dna') +
      facet_wrap(~seq_group, ncol=1) + xlab('Position') + ylab('gkmexplain importance score')
    pdf(gsub(':', '_', glue('{SNP}.{CLUSTER}.pdf')), height=8, width=10)
    print(p)
    dev.off()
  }
}
