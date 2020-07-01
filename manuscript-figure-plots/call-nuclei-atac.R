#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)
library(viridis)
library(ggpointdensity)

parse_nucleus <- function(x) {
  RE <- '(.*)-(.*)'
  lib <- gsub(RE, '\\1', x)
  barcode <- gsub(RE, '\\2', x)
  return(list('library'=lib, 'barcode'=barcode))
}

args <- commandArgs(T)
METRICS <- args[1]
INITIAL_THRESHOLDS <- args[2]
FINAL_THRESHOLDS <- args[3]
LIBRARY_LABELS <- args[4]

library_labels <- read.table(LIBRARY_LABELS, head=T, sep='\t', stringsAsFactors = F)
library_labels <- library_labels[library_labels$modality=='ATAC' & library_labels$fans_status=='no FANS',]

thresholds <- read.table(INITIAL_THRESHOLDS, head=T, as.is=T, sep='\t')
final_thresholds <- read.table(FINAL_THRESHOLDS, head=T, as.is=T, sep='\t')
thresholds <- left_join(thresholds %>% dplyr::select(-max_hqaa), final_thresholds %>% dplyr::rename(max_hqaa=max_hqaa_threshold))

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
assign_species$drop <- apply(assign_species[,c('library', 'barcode', 'best', 'worst', 'assignment')], 1, function(x){
  if (x[5]=='none') {
    return(paste(paste(x[1], genomes, x[2], sep='-'), collapse=','))
  } else {
    return(paste(paste(x[1], genomes[genomes!=x[3]], x[2], sep='-'), collapse=','))
  }
})
drop <- unlist(strsplit(assign_species$drop, ','))

ataqv <- ataqv[!ataqv$nucleus %in% drop,]


# Plot for all libraries we used downstream...
used_downstream <- ataqv[ataqv$library %in% library_labels$library,]
used_downstream <- left_join(used_downstream, library_labels[,c('library', 'name')])
used_downstream_thresholds <- unique(used_downstream[,c('library', 'min_hqaa', 'max_hqaa', 'max_max_fraction_reads_from_single_autosome', 'min_tss_enrichment', 'name')])

p <- ggplot(used_downstream) + geom_point(aes(x = hqaa, y = tss_enrichment), alpha=0.02, stroke=0) +
#p <- ggplot(used_downstream) + geom_pointdensity(aes(x = hqaa, y = tss_enrichment), size=1) +
  facet_wrap(~name) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  geom_vline(aes(xintercept = min_hqaa), col='red', linetype='dashed', data = used_downstream_thresholds) +
  geom_vline(aes(xintercept = max_hqaa), col='red', linetype='dashed', data = used_downstream_thresholds) +
  geom_hline(aes(yintercept = min_tss_enrichment), col='red', linetype='dashed', data=used_downstream_thresholds) +
  xlab('Final high-quality autosomal reads') +
  ylab('TSS enrichment') +
  theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
  scale_color_viridis(trans='log') +
  guides(color=guide_colorbar(title='log(n_neighbors)'))
png('hqaa-vs-tss-enrichment-used-downstream.png', height = 4, width = 7, res=300, units='in')
p
dev.off()

p <- ggplot(used_downstream) + geom_point(aes(x = hqaa, y = max_fraction_reads_from_single_autosome), alpha=0.05, stroke=0) +
  facet_wrap(~name) +
  scale_x_log10() +
  theme_bw() +
  geom_vline(aes(xintercept = min_hqaa), col='red', linetype='dashed', data = used_downstream_thresholds) +
  geom_vline(aes(xintercept = max_hqaa), col='red', linetype='dashed', data = used_downstream_thresholds) +
  geom_hline(aes(yintercept = max_max_fraction_reads_from_single_autosome), col='red', linetype='dashed', data=used_downstream_thresholds) +
  xlab('Final high-quality autosomal reads') +
  ylab('Max. fraction reads from single autosome') +
  ylim(c(0, 0.5)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
  scale_color_viridis(trans='log') +
  guides(color=guide_colorbar(title='log(n_neighbors)'))
png('hqaa-vs-max-fraction-reads-from-single-autosome-used-downstream.png', height = 4, width = 7, res=300, units='in')
p
dev.off()


# First, do the plots for the FANS vs no FANS
fans_vs_no_fans <- ataqv[ataqv$library %in% c('133151-hg19', '133152-hg19', '133153-hg19', '133154-hg19'),] %>%
  #dplyr::mutate(fans_status=ifelse(library %in% c('133151-hg19', '133153-hg19'), 'Crude', 'FANS'),
  dplyr::mutate(fans_status=ifelse(library %in% c('133151-hg19', '133153-hg19'), 'no FANS', 'FANS'),
                replicate=ifelse(library %in% c('133151-hg19', '133152-hg19'), 'rep. 1', 'rep. 2'))
fans_vs_no_fans_thresholds <- unique(fans_vs_no_fans[,c('library', 'min_hqaa', 'max_hqaa', 'max_max_fraction_reads_from_single_autosome', 'min_tss_enrichment', 'fans_status', 'replicate')])

p <- ggplot(fans_vs_no_fans) + geom_point(aes(x = hqaa, y = tss_enrichment), alpha=0.05, stroke=0) +
  facet_grid(replicate~fans_status) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  geom_vline(aes(xintercept = min_hqaa), col='red', linetype='dashed', data = fans_vs_no_fans_thresholds) +
  geom_vline(aes(xintercept = max_hqaa), col='red', linetype='dashed', data = fans_vs_no_fans_thresholds) +
  geom_hline(aes(yintercept = min_tss_enrichment), col='red', linetype='dashed') +
  xlab('Final high-quality autosomal reads') +
  ylab('TSS enrichment')
png('hqaa-vs-tss-enrichment-fans-vs-no-fans.png', height = 5, width = 5, res=300, units='in')
p
dev.off()

p <- ggplot(fans_vs_no_fans) + geom_point(aes(x = hqaa, y = max_fraction_reads_from_single_autosome), alpha=0.05, stroke=0) +
  facet_grid(replicate~fans_status) +
  scale_x_log10() +
  theme_bw() +
  geom_vline(aes(xintercept = min_hqaa), col='red', linetype='dashed', data = fans_vs_no_fans_thresholds) +
  geom_vline(aes(xintercept = max_hqaa), col='red', linetype='dashed', data = fans_vs_no_fans_thresholds) +
  geom_hline(aes(yintercept = max_max_fraction_reads_from_single_autosome), col='red', linetype='dashed') +
  xlab('Final high-quality autosomal reads') +
  ylab('Max. fraction reads from single autosome') +
  ylim(c(0, 0.5))
png('hqaa-vs-max-fraction-reads-from-single-autosome-fans-vs-no-fans.png', height = 5, width = 5, res=300, units='in')
p
dev.off()



# then, do the plots for the 20k vs 40k
loading <- ataqv[ataqv$library %in% c('63_20-hg19', '63_40-hg19'),] %>%
  dplyr::mutate(nuclei_loaded=paste0(gsub('63_(\\d+)-hg19', '\\1', library), 'k nuclei'))
loading_thresholds <- unique(loading[,c('library', 'min_hqaa', 'max_hqaa', 'max_max_fraction_reads_from_single_autosome', 'min_tss_enrichment', 'nuclei_loaded')])

p <- ggplot(loading) + geom_point(aes(x = hqaa, y = tss_enrichment), alpha=0.05, stroke=0) +
  facet_wrap(~nuclei_loaded) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  geom_vline(aes(xintercept = min_hqaa), col='red', linetype='dashed', data = loading_thresholds) +
  geom_vline(aes(xintercept = max_hqaa), col='red', linetype='dashed', data = loading_thresholds) +
  geom_hline(aes(yintercept = min_tss_enrichment), col='red', linetype='dashed') +
  xlab('Final high-quality autosomal reads') +
  ylab('TSS enrichment')
png('hqaa-vs-tss-enrichment-20k-vs-40k.png', height = 3, width = 5, res=300, units='in')
p
dev.off()

p <- ggplot(loading) + geom_point(aes(x = hqaa, y = max_fraction_reads_from_single_autosome), alpha=0.05, stroke=0) +
  facet_wrap(~nuclei_loaded) +
  scale_x_log10() +
  theme_bw() +
  geom_vline(aes(xintercept = min_hqaa), col='red', linetype='dashed', data = loading_thresholds) +
  geom_vline(aes(xintercept = max_hqaa), col='red', linetype='dashed', data = loading_thresholds) +
  geom_hline(aes(yintercept = max_max_fraction_reads_from_single_autosome), col='red', linetype='dashed') +
  xlab('Final high-quality autosomal reads') +
  ylab('Max. fraction reads from single autosome') +
  ylim(c(0, 0.5))
png('hqaa-vs-max-fraction-reads-from-single-autosome-20k-vs-40k.png', height = 4, width = 6, res=300, units='in')
p
dev.off()



# then, do the plots for the remaining rat and human library
feb <- ataqv[ataqv$library %in% c('125589-hg19', '125589-rn6'),] %>%
  dplyr::mutate(species=ifelse(gsub('125589-', '', library)=='hg19', 'human', 'rat'))
feb_thresholds <- unique(feb[,c('library', 'min_hqaa', 'max_hqaa', 'max_max_fraction_reads_from_single_autosome', 'min_tss_enrichment', 'species')])

p <- ggplot(feb) + geom_point(aes(x = hqaa, y = tss_enrichment), alpha=0.05, stroke=0) +
  facet_wrap(~species) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  geom_vline(aes(xintercept = min_hqaa), col='red', linetype='dashed', data = feb_thresholds) +
  geom_vline(aes(xintercept = max_hqaa), col='red', linetype='dashed', data = feb_thresholds) +
  geom_hline(aes(yintercept = min_tss_enrichment), col='red', linetype='dashed') +
  xlab('Final high-quality autosomal reads') +
  ylab('TSS enrichment')
png('hqaa-vs-tss-enrichment-feb.png', height = 3, width = 5, res=300, units='in')
p
dev.off()

p <- ggplot(feb) + geom_point(aes(x = hqaa, y = max_fraction_reads_from_single_autosome), alpha=0.05, stroke=0) +
  facet_wrap(~species) +
  scale_x_log10() +
  theme_bw() +
  geom_vline(aes(xintercept = min_hqaa), col='red', linetype='dashed', data = feb_thresholds) +
  geom_vline(aes(xintercept = max_hqaa), col='red', linetype='dashed', data = feb_thresholds) +
  geom_hline(aes(yintercept = max_max_fraction_reads_from_single_autosome), col='red', linetype='dashed') +
  xlab('Final high-quality autosomal reads') +
  ylab('Max. fraction reads from single autosome') +
  ylim(c(0, 0.5))
png('hqaa-vs-max-fraction-reads-from-single-autosome-feb.png', height = 3, width = 5, res=300, units='in')
p
dev.off()



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
pdf('hqaa-vs-tss-enrichment.pdf', height = 10, width = 10)
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
pdf('hqaa-vs-max-fraction-reads-from-single-autosome.pdf', height = 10, width = 10)
p
dev.off()


survivors <- ataqv %>%
  dplyr::filter(hqaa>=min_hqaa) %>%
  dplyr::filter(hqaa<=max_hqaa) %>%
  dplyr::filter(tss_enrichment>=min_tss_enrichment) %>%
  dplyr::filter(max_fraction_reads_from_single_autosome<=max_max_fraction_reads_from_single_autosome)
survivors$barcode <- parse_nucleus(survivors$nucleus)$barcode

write.table(survivors %>% dplyr::select(library, barcode), 'atac-nuclei.txt', append = F, quote = F, sep='\t', row.names = F, col.names = F)

# also write out the table of thresholds used...
supplemental_table_of_thresholds_used <- left_join(library_labels, thresholds) %>%
  dplyr::select(library, experiment, min_hqaa, max_hqaa, min_tss_enrichment, max_max_fraction_reads_from_single_autosome)
supplemental_table_of_thresholds_used$max_hqaa <- floor(supplemental_table_of_thresholds_used$max_hqaa)
supplemental_table_of_thresholds_used$library[supplemental_table_of_thresholds_used$library=='125589-hg19'] <- '122589 (human nuclei)'
supplemental_table_of_thresholds_used$library[supplemental_table_of_thresholds_used$library=='125589-rn6'] <- '122589 (rat nuclei)'
colnames(supplemental_table_of_thresholds_used) <- c('library', 'experiment', 'min. reads after filtering', 'max. reads after filtering', 'min. TSS enrichment', 'max(max fraction reads from single autosome)')
write.table(supplemental_table_of_thresholds_used, file = 'atac-qc-thresholds.csv', append=F, quote = F, sep = ',', row.names = F, col.names = T)
