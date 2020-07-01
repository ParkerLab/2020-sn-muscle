#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)

args <- commandArgs(T)
BINARY_DHS_MATRIX <- args[1]
METADATA <- args[2]

load(BINARY_DHS_MATRIX) # object: dat_bin

sample_info <- read.table(METADATA, head=T, as.is=T, sep='\t', quote = '', comment.char = '')
sample_info$library_name <- paste(sample_info$Biosample.name, sample_info$Altius.Biosample.ID, sep='.')
stopifnot(all(sample_info$library_name==colnames(dat_bin)))

sample_info[sample_info$Organ=='' | is.na(sample_info$Organ),]
sample_info <- sample_info[sample_info$Organ!='',]

# are there any tissues where the only samples are cancer?
unique(sample_info[,c('Biosample.type', 'Biological.state')])

tmp <- sample_info %>%
  dplyr::group_by(Organ, Biosample.type) %>%
  dplyr::summarize(count=n(),
                   have_samples=1)

p <- ggplot(tmp) +
  geom_tile(aes(x=Organ, y=Biosample.type, fill=have_samples)) +
  theme_bw() +
  scale_fill_gradient(high='red', low='white') +
  coord_flip()
#p

ALL_TISSUES <- unique(sample_info$Organ)
HAVE_PRIMARY_SAMPLES <- unique(sample_info$Organ[sample_info$Biosample.type=='Primary'])
MISSING_PRIMARY_SAMPLES <- ALL_TISSUES[!ALL_TISSUES %in% HAVE_PRIMARY_SAMPLES]
print('Missing primary samples:')
print(paste(MISSING_PRIMARY_SAMPLES, collapse=', '))

# What is the difference between Umbilical and Umbilical Cord?
# sample_info[grep('Umbilical', sample_info$Organ),]

# Drop some, especially those that are similar to our cell types/adipose/beta cells for which we have ATAC
DROP_TISSUES <- c('Blood', 'Pulmonary Artery', 'Umbilical', 'Vascular', 'Muscle', 'Stroma', 'Pancreas', 'Tongue', 'Adipose')
KEEP <- HAVE_PRIMARY_SAMPLES[!HAVE_PRIMARY_SAMPLES %in% DROP_TISSUES]

sample_info <- sample_info[sample_info$Biosample.type=='Primary' & sample_info$Organ %in% KEEP,]
table(sample_info$Organ)

for(tissue in KEEP) {
  samples <- sample_info$library_name[sample_info$Organ==tissue]
  samples_dhs <- dat_bin[,samples]
  if (length(samples) > 1) {
    number_samples_with_peak <- apply(samples_dhs, 1, sum)
    #tissue_peaks <- data.frame(peak=names(number_samples_with_peak[number_samples_with_peak>=floor((length(samples)/2)+1)]))
    tissue_peaks <- data.frame(peak=names(number_samples_with_peak[number_samples_with_peak>=length(samples)/2]))
    write.table(tissue_peaks, file=gsub(' ', '_', glue('{tissue}.hg38-peaks.txt')), append = F, quote = F, sep = '\t', row.names = F, col.names = F)
  } else {
    tissue_peaks <- data.frame(peak=names(samples_dhs[samples_dhs==T]))
    write.table(tissue_peaks, file=gsub(' ', '_', glue('{tissue}.hg38-peaks.txt')), append = F, quote = F, sep = '\t', row.names = F, col.names = F)
  }
}
