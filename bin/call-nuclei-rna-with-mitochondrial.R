#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(glue)
library(tidyr)

args <- commandArgs(T)
THRESHOLDS <- args[1]
QC <- args[2:length(args)]

# for testing
# QC <- list.files('/lab/work/porchard/sn-muscle-project/work/downstream/results/nucleus-qc', pattern='hg19', full.names = T)


qc <- bind_rows(lapply(QC, function(f){
  lib <- gsub('.txt', '', basename(f))
  tmp <- read.table(f, head = T, sep = '\t', as.is = T)
  tmp$library <- lib
  return(tmp)
}))

#thresholds <- data.frame(
#  library=c('133155-hg19', '133156-hg19', '133157-hg19', '133158-hg19', '63_20_rna-hg19', '63_40_rna-hg19'),
#  min_hqaa=c(5000, 3000, 5000, 3000, 1000, 1000),
#  max_hqaa=c(35000, 25000, 35000, 25000, 9000, 9000),
#  max_mitochondrial=c(0.003, 0.007, 0.003, 0.007, 0.005, 0.005)
#)
thresholds <- read.table(THRESHOLDS, stringsAsFactors = F, sep='\t', header = T)

p <- ggplot(qc) +
  geom_point(aes(x = number_umis, y=fraction_mitochondrial+0.001), alpha=0.1) +
  geom_vline(aes(xintercept = min_hqaa), data = thresholds, color='red', linetype='dashed') +
  geom_vline(aes(xintercept = max_hqaa), data = thresholds, color='red', linetype='dashed') +
  geom_hline(aes(yintercept = max_mitochondrial), data = thresholds, color='red', linetype='dashed') +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~library) +
  theme_bw() +
  xlab('Number of UMIs') +
  ylab('Fraction mitochondrial + 0.001')
png('hqaa-vs-mitochondrial.png', height=6, width=12, units='in', res=300)
p
dev.off()

pass_qc <- qc %>%
  left_join(thresholds) %>%
  dplyr::filter(number_umis>=min_hqaa) %>%
  dplyr::filter(number_umis<=max_hqaa) %>%
  dplyr::filter(fraction_mitochondrial+0.001<=max_mitochondrial) %>%
  dplyr::select(library, barcode)

hqaa <- qc %>%
  left_join(thresholds) %>%
  dplyr::filter(number_umis>=min_hqaa) %>%
  dplyr::filter(number_umis<=max_hqaa) %>%
  dplyr::filter(fraction_mitochondrial+0.001<=max_mitochondrial) %>%
  dplyr::select(library, barcode, number_umis)

write.table(pass_qc, file='rna-nuclei.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = F)
write.table(hqaa, file='hqaa.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = F)
