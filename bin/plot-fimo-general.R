#!/usr/bin/env Rscript
library(ggplot2)
library(tidyr)
library(dplyr)
library(glue)
library(ggseqlogo)

args <- commandArgs(T)
MOTIF_DIR <- args[1]
GKMEXPLAIN_DIR <- args[2]
SEQUENCE_NAME <- args[3] # e.g. 2:23934533:G:A
PLAIN_MOTIF_DIR <- args[4]

EXPLAINED_FILES <- list.files(GKMEXPLAIN_DIR, full.names=T, pattern='reformatted.txt')
MOTIF_FILES <- list.files(MOTIF_DIR, pattern = '.fimo.txt', full.names=T)
SEQUENCE_NAME <- gsub('chr', '', SEQUENCE_NAME)

load_fimo_file <- function(f) {
  FILENAME_RE <- '(.*)\\.fimo\\.txt$'
  tmp <- read.table(f, head=T, stringsAsFactors = F, sep='\t')
  tmp <- tmp[tmp$motif_id != 'motif_id',c('motif_id', 'sequence_name', 'start', 'stop', 'strand', 'score', 'p.value', 'q.value', 'matched_sequence')]
  tmp$start <- as.numeric(tmp$start)
  tmp$stop <- as.numeric(tmp$stop)
  tmp$score <- as.numeric(tmp$score)
  tmp$p.value <- as.numeric(tmp$p.value)
  tmp$ref_or_alt <- gsub(FILENAME_RE, '\\1', basename(f))
  return(tmp)
}

fimo <- bind_rows(lapply(MOTIF_FILES, load_fimo_file)) %>% dplyr::select(motif_id, sequence_name, start, stop, strand, p.value, ref_or_alt) %>%
  tidyr::spread(key=ref_or_alt, value=p.value, fill=1) %>%
  dplyr::filter(alt!=ref)

# for each of the motifs, plot their pwms against the motif discovered by gkmexplain
# load in the gkmexplain results

read_explain_file <- function(f) {
  FILENAME_RE <- '(cluster_.*)\\.(ref|alt).explained.reformatted.txt'
  tmp <- read.table(f, head=F, sep='\t', stringsAsFactors = F)
  colnames(tmp) <- c('snp', 'pos', 'nuc', 'score')
  tmp$cluster <- gsub(FILENAME_RE, '\\1', basename(f))
  tmp$ref_alt <- gsub(FILENAME_RE, '\\2', basename(f))
  return(tmp)
}


explained <- bind_rows(lapply(EXPLAINED_FILES, read_explain_file))
explained$snp <- gsub('chr', '', explained$snp)

fimo$sequence_name <- gsub('chr', '', fimo$sequence_name) 
fimo <- fimo[fimo$sequence_name %in% explained$snp,]
fimo <- fimo[fimo$sequence_name==SEQUENCE_NAME,]

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

#for(SNP in unique(fimo$sequence_name)) {
  for(CLUSTER in unique(explained$cluster)) {
    for (i in 1:nrow(fimo)) {
      MOTIF <- fimo$motif_id[i]
      FIMO_MOTIF_START <- fimo$start[i]
      STRAND <- fimo$strand[i]
      SNP <- fimo$sequence_name[i]
      
      # position 26 in gkmexplain should be position 51 in fimo...
      motif_file <- gsub('::', '_', glue('{PLAIN_MOTIF_DIR}/{MOTIF}.txt'))
      if (!file.exists(motif_file)){
        next
      }
      motif <- as.matrix(read.table(motif_file, head=F, as.is=F, sep='\t'))
      colnames(motif) <- seq(0, ncol(motif)-1) + FIMO_MOTIF_START - 51 + 26
      rownames(motif) <- c('A', 'C', 'G', 'T')
      if(max(as.numeric(colnames(motif))) > 51) {
        next
      }
      if(min(as.numeric(colnames(motif))) < 1) {
        next
      }
      if (STRAND=='-') {
        motif <- motif[,rev(colnames(motif))]
        motif <- motif[rev(rownames(motif)),]
        colnames(motif) <- rev(colnames(motif))
        rownames(motif) <- c('A', 'C', 'G', 'T')
      }
      motif <- cbind(matrix(0, ncol=as.numeric(colnames(motif)[1])-1, nrow=4), motif, matrix(0, ncol=51-as.numeric(colnames(motif)[ncol(motif)]), nrow=4))
      colnames(motif) <- 1:51
      
      gkmexplain_ref <- to_matrix(CLUSTER, SNP, 'ref', positions=50)
      gkmexplain_alt <- to_matrix(CLUSTER, SNP, 'alt', positions=50)
      ref_minus_alt <- gkmexplain_ref - gkmexplain_alt
      p <- ggseqlogo(list('With ref. allele'=gkmexplain_ref, 'With alt. allele'=gkmexplain_alt, 'Ref - alt'=ref_minus_alt), method='custom', seq_type='dna') +
        facet_wrap(~seq_group, ncol=1) + xlab('') + ylab('gkmexplain importance score') +
        theme(axis.text.x = element_blank())
      p2 <- ggseqlogo(motif, method='bits', seq_type='dna') +
        xlab('') +
        ggtitle(MOTIF) +
        theme(axis.text.x = element_blank())
      STRAND <- ifelse(STRAND=='+', 'fwd', 'rev')
      pdf(glue('{SNP}___{CLUSTER}___{MOTIF}_{STRAND}_{FIMO_MOTIF_START}.pdf'))
      print(cowplot::plot_grid(plotlist=list(p, p2), ncol = 1, rel_heights = c(3, 1), align = 'v'))
      dev.off()
    }
  }
#}
