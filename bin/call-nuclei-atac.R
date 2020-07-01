#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)

parse_nucleus <- function(x) {
  RE <- '(.*)-(.*)'
  lib <- gsub(RE, '\\1', x)
  barcode <- gsub(RE, '\\2', x)
  return(list('library'=lib, 'barcode'=barcode))
}

args <- commandArgs(T)
METRICS <- args[1]
THRESHOLDS <- args[2]

#thresholds <- data.frame(
#  library=c('125589-hg19', '125589-rn6', '133151-hg19', '133152-hg19', '133153-hg19', '133154-hg19', '63_20-hg19', '63_40-hg19'),
#  min_hqaa=c(50000, 50000, 70000, 20000, 50000, 20000, 25000, 15000),
#  max_hqaa=c(170000, 170000, 200000, 130000, 200000, 130000, 120000, 80000),
#  max_max_fraction_reads_from_single_autosome=c(0.13, 0.15, rep(0.13, 6)),
#  min_tss_enrichment=c(3, 3, 4.5, 3, 4.5, 3, 4.5, 4.5))
thresholds <- read.table(THRESHOLDS, stringsAsFactors = F, sep='\t', header = T)


ataqv <- read.table(METRICS, head = F, as.is = T, sep = '\t', col.names = c('nucleus', 'metric', 'value'), colClasses = c('character', 'character', 'character'))
ataqv$library <- parse_nucleus(ataqv$nucleus)$library
ataqv$value[ataqv$value=='None' | ataqv$value=='NA' | is.na(ataqv$value)] <- NA
ataqv$value <- as.numeric(ataqv$value)
ataqv <- ataqv[ataqv$metric %in% c('tss_enrichment', 'hqaa', 'max_fraction_reads_from_single_autosome'),]
ataqv <- ataqv %>%
  tidyr::spread(key = metric, value = value) %>%
  left_join(thresholds)

# assign species where appropriate...
assign_species <- ataqv
assign_species$genome <- gsub('(.*)-(.*)', '\\2', assign_species$library)
assign_species$library <- gsub('(.*)-(.*)', '\\1', assign_species$library)
number_species_per_library <- assign_species %>%
  dplyr::group_by(library) %>%
  dplyr::summarize(number_species=n_distinct(genome))
libraries_with_multiple_species <- number_species_per_library$library[number_species_per_library$number_species>1]
assign_species <- assign_species[assign_species$library %in% libraries_with_multiple_species,]
assign_species$barcode <- gsub('.*-', '', assign_species$nucleus)
assign_species <- assign_species %>%
  dplyr::select(library, genome, barcode, hqaa) %>%
  tidyr::spread(key=genome, value=hqaa)
genomes <- colnames(assign_species)[3:ncol(assign_species)]
assign_species$ratio <- apply(assign_species[,genomes], 1, function(x){max(x/sum(x))})
assign_species$best <- apply(assign_species[,genomes], 1, function(x){genomes[x==max(x)][1]})
assign_species$worst <- apply(assign_species[,genomes], 1, function(x){genomes[x==min(x)][1]})
assign_species$assignment <- 'none'
assign_species$assignment[assign_species$ratio>=0.87] <- assign_species$best[assign_species$ratio>=0.87]
p <- ggplot(assign_species) +
  geom_point(aes(x=hg19+rn6, y=ratio), alpha=0.2, stroke=0) +
  scale_x_log10() +
  ylab('Max(hg19/(rn6+hg19), rn6/(rn6+hg19))') +
  xlab('hg19+rn6') +
  theme_bw() +
  geom_hline(yintercept = 0.87, linetype='dashed', color='red')
png('hg19-rn6-ratio-threshold.png', height=5, width=6)
p
dev.off()
assign_species$drop <- apply(assign_species[,c('library', 'barcode', 'best', 'worst', 'assignment')], 1, function(x){
  if (x[5]=='none') {
    return(paste(paste(x[1], genomes, x[2], sep='-'), collapse=','))
  } else {
    return(paste(paste(x[1], genomes[genomes!=x[3]], x[2], sep='-'), collapse=','))
  }
})
drop <- unlist(strsplit(assign_species$drop, ','))

ataqv <- ataqv[!ataqv$nucleus %in% drop,]


# will do this on a per-library basis...
p <- ggplot(ataqv) + geom_point(aes(x = hqaa, y = tss_enrichment), alpha=0.05, stroke=0) +
  facet_wrap(~library) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  geom_vline(aes(xintercept = min_hqaa), col='red', linetype='dashed', data = thresholds) +
  geom_vline(aes(xintercept = max_hqaa), col='red', linetype='dashed', data = thresholds) +
  geom_hline(aes(yintercept = min_tss_enrichment), col='red', linetype='dashed') +
  xlab('Final high-quality autosomal reads') +
  ylab('TSS enrichment')
png('hqaa-vs-tss-enrichment.png', height = 10, width = 10, units='in', res=300)
p
dev.off()

p <- ggplot(ataqv) + geom_point(aes(x = hqaa, y = max_fraction_reads_from_single_autosome), alpha=0.05, stroke=0) +
  facet_wrap(~library) +
  scale_x_log10() +
  theme_bw() +
  geom_vline(aes(xintercept = min_hqaa), col='red', linetype='dashed', data = thresholds) +
  geom_vline(aes(xintercept = max_hqaa), col='red', linetype='dashed', data = thresholds) +
  geom_hline(aes(yintercept = max_max_fraction_reads_from_single_autosome), col='red', linetype='dashed') +
  xlab('Final high-quality autosomal reads') +
  ylab('Max. fraction reads from single autosome')
png('hqaa-vs-max-fraction-reads-from-single-autosome.png', height = 10, width = 10, res=300, units='in')
p
dev.off()


survivors <- ataqv %>%
  dplyr::filter(hqaa>=min_hqaa) %>%
  dplyr::filter(hqaa<=max_hqaa) %>%
  dplyr::filter(tss_enrichment>=min_tss_enrichment) %>%
  dplyr::filter(max_fraction_reads_from_single_autosome<=max_max_fraction_reads_from_single_autosome)
survivors$barcode <- parse_nucleus(survivors$nucleus)$barcode

write.table(survivors %>% dplyr::select(library, barcode), 'atac-nuclei.txt', append = F, quote = F, sep='\t', row.names = F, col.names = F)
