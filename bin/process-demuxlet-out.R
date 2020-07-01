#!/usr/bin/env Rscript
library(ggplot2)
library(tidyr)
library(dplyr)
library(glue)

# for testing
# BEST_FILES <- list.files('/lab/work/porchard/sn-muscle-project/work/downstream-final/results/demuxlet', full.names = T, pattern='.best')

args <- commandArgs(T)
BEST_FILES <- args

load_best_file <- function(f) {
  tmp <- read.table(f, head=T, sep='\t')
  tmp$assignment <- 'doublet'
  tmp$assignment[grep('SNG-KSM1', tmp$BEST)] <- 'KSM1'
  tmp$assignment[grep('SNG-KSM2', tmp$BEST)] <- 'KSM2'
  tmp$library <- gsub('.best', '', basename(f))
  return(tmp)
}

best <- bind_rows(lapply(BEST_FILES, load_best_file))
best <- best[,c('library', 'BARCODE', 'assignment')]

counts <- best %>%
  dplyr::group_by(library, assignment) %>%
  dplyr::summarize(count=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(library) %>%
  dplyr::mutate(fraction=count/sum(count),
                perc=round(fraction*100, 1)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(label=glue('{count} ({perc}%)'))

counts$library <- gsub('63_20-hg19', '20k input', counts$library)
counts$library <- gsub('63_40-hg19', '40k input', counts$library)

p <- ggplot(counts) +
  geom_bar(aes(x=assignment, y=fraction, fill=library), position=position_dodge(width=1), stat='identity') +
  geom_text(aes(x=assignment, y=fraction+0.01, label=label, fill=library), position=position_dodge(width=1)) +
  theme_bw() +
  ylab('Fraction of pass-QC nuclei') +
  xlab('Sample or doublet') +
  guides(fill=guide_legend(title=''))
pdf('snATAC-doublets.pdf', height=5, width=3)
print(p)
dev.off()

write.table(best, file='atac-demuxlet-assignments.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = F)
