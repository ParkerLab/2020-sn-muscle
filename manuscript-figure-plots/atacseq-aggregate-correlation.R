#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(GGally)

args <- commandArgs(T)
LIBRARY_LABELS <- args[1]
COUNT_FILES <- args[2:length(args)]
counts <- bind_rows(lapply(COUNT_FILES, function(f){tmp <- read.table(f, head = T, as.is = T, sep = '\t', col.names = c('chrom', 'start', 'end', 'count')); tmp$f <- f; return(tmp)}))
counts$library <- gsub('.counts.bed', '', basename(counts$f))
counts <- counts %>%
  dplyr::select(-f) %>%
  tidyr::spread(key=library, value=count)
lengths <- counts$end - counts$start
counts <- counts %>%
  dplyr::select(-chrom, -start, -end)

# convert to rpkm...
counts_to_rpkm <- function(counts_vector, lengths) {
  # lengths indicates the lengths (in bp) of each of the corresponding regions
  stopifnot(length(counts_vector) == length(lengths))
  
  counts_sum <- sum(counts_vector)
  
  # to reduce the probability of integer overflow,
  # enforce a certain order of operations
  rpkm <- counts_vector * (((10^9) / (lengths)) / counts_sum)
  return(rpkm)
}

for(i in colnames(counts)) {
  counts[,i] <- counts_to_rpkm(counts[,i], lengths)
}

colnames(counts) <- gsub('.pruned', '', colnames(counts))

# plot correlation for all libraries present
library_labels <- read.table(LIBRARY_LABELS, head = T, as.is = T, sep = '\t')
fans_info <- library_labels[!is.na(library_labels$experiment) & library_labels$experiment=='FANS vs. no FANS',]
#fans_info$fans_status <- ifelse(grepl('no FANS', fans_info$name), 'Crude', 'FANS')
fans_info$fans_status <- ifelse(grepl('no FANS', fans_info$name), 'no FANS', 'FANS')
fans_info$label <- paste(fans_info$fans_status, fans_info$replicate)
loading_info <- library_labels[!is.na(library_labels$experiment) & library_labels$experiment=='20k vs. 40k',]
loading_info$loading_concentration <- ifelse(grepl('20', loading_info$name), '20k nuclei', '40k nuclei')
loading_info$label <- loading_info$loading_concentration
bulk_info <- data.frame(library=c('320-NM-1', '320-NM-2', '320-NM-3', '320-NM-4'), label=c('HSM1 bulk rep. 1', 'HSM1 bulk rep. 2', 'HSM2 bulk rep. 1', 'HSM2 bulk rep. 2'))
all_info <- bind_rows(fans_info[,c('library', 'label')], bulk_info, loading_info[,c('library', 'label')])

all <- counts[,colnames(counts) %in% all_info$library]
colnames(all) <- sapply(colnames(all), function(x){all_info$label[all_info$library==x]})
all <- log2(all+0.01)

p <- ggpairs(all, lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw() +
  xlab('log2(RPKM + 0.01)') +
  ylab('log2(RPKM + 0.01)')
png('all-atac-correlation.png', width = 10, height = 10, units='in', res=300)
p
dev.off()
