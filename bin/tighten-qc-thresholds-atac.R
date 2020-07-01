#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)
library(viridis)

# tighten up ATAC QC thresholds given the output of demuxlet and the HQAA / TSS enrichment statistics...
#DEMUXLET <- '/lab/work/porchard/sn-muscle-project/work/downstream-final/results/demuxlet/atac-demuxlet-assignments.txt'
#ATAC_NUCLEI <- '/lab/work/porchard/sn-muscle-project/work/downstream-final/results/nucleus-qc/atac-nuclei.txt'
#ATAQV_METRICS <- '/lab/work/porchard/sn-muscle-project/work/downstream-final/results/nucleus-qc/metrics.txt'
args <- commandArgs(T)
DEMUXLET <- args[1]
ATAC_NUCLEI <- args[2]
ATAQV_METRICS <- args[3]

MIXED_SPECIES_LIBRARIES <- c('125589-hg19', '125589-rn6') # half of the doublets have already been filtered out of these, so need to take that into account
nuclei <- read.table(ATAC_NUCLEI, stringsAsFactors = F, sep='\t', col.names=c('library', 'barcode'))

nuclei_assignments <- bind_rows(lapply(c(DEMUXLET), function(f){
  tmp <- read.table(f, head=F, sep='\t', col.names = c('library', 'barcode', 'assignment'), colClasses = c('character'))
  return(tmp)
}))

ataqv <- read.table(ATAQV_METRICS, sep='\t', stringsAsFactors = F, col.names = c('nucleus', 'metric', 'value')) %>%
  dplyr::filter(metric %in% c('hqaa', 'tss_enrichment', 'percent_duplicate')) %>%
  dplyr::mutate(value=as.numeric(value),
                library=gsub('(.*)-(.*)', '\\1', nucleus),
                barcode=gsub('(.*)-(.*)', '\\2', nucleus)) %>%
  dplyr::select(-nucleus) %>%
  tidyr::spread(key=metric, value=value)
ataqv <- ataqv %>%
  dplyr::semi_join(nuclei)

atac <- left_join(nuclei_assignments, ataqv) %>%
  dplyr::filter(!is.na(hqaa))

# plot out log10(HQAA) vs local probability of doublet (in 100-nucleus bins?)
binned <- atac[order(atac$library, atac$hqaa),] %>%
  dplyr::group_by(library) %>%
  dplyr::mutate(bin=ntile(n=n()/100)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(library, bin) %>%
  dplyr::summarize(median_hqaa=median(hqaa),
                   doublets=sum(assignment=='doublet'),
                   singlets=sum(assignment!='doublet'),
                   prob_doublet=doublets/(doublets+singlets))


# the doublet rate is 1/2 the actual value
# if the cutoff point should be where the probability that you're a doublet surpasses 0.5,
# then we want to find the first bin where at least 1/4 of the barcodes were assigned as doublets
thresholds <- binned %>%
  dplyr::filter((prob_doublet*2)>=0.5) %>%
  dplyr::group_by(library) %>%
  dplyr::summarize(max_hqaa_threshold=min(median_hqaa))

p <- ggplot(binned) +
  geom_point(aes(x=median_hqaa, y=prob_doublet, color=library)) +
  geom_line(aes(x=median_hqaa, y=prob_doublet, color=library)) +
  geom_vline(aes(xintercept=max_hqaa_threshold, color=library), data=thresholds, linetype='dashed') +
  theme_bw()
pdf('hqaa-vs-doublet-assignments.pdf', height=4, width=7)
p
dev.off()

# use the quantiles to estimate the proper thresholds...?
quantiles_by_all <- left_join(ataqv, thresholds) %>%
  dplyr::filter(!is.na(max_hqaa_threshold)) %>%
  dplyr::mutate(remove=hqaa>max_hqaa_threshold) %>%
  dplyr::group_by(library) %>%
  dplyr::summarize(fraction_removed=mean(remove))

quantiles_by_singlets <- left_join(nuclei_assignments, ataqv) %>%
  left_join(thresholds) %>%
  dplyr::filter(!is.na(max_hqaa_threshold)) %>%
  dplyr::filter(assignment!='doublet') %>%
  dplyr::mutate(remove=hqaa>max_hqaa_threshold) %>%
  dplyr::group_by(library) %>%
  dplyr::summarize(fraction_removed=mean(remove))


# the other libraries are ~15k nuclei, which should resemble the 20k nuclei doublet rate
# so take that as the fraction to remove for those libraries...
quantile_for_unmixed <- 1-quantiles_by_all$fraction_removed[quantiles_by_all$library=='63_20-hg19']
thresholds_for_others <- ataqv[!ataqv$library %in% c(quantiles_by_all$library, MIXED_SPECIES_LIBRARIES),] %>%
  dplyr::group_by(library) %>%
  dplyr::summarize(max_hqaa_threshold=quantile(hqaa, c(quantile_for_unmixed)))
thresholds <- bind_rows(thresholds, thresholds_for_others)

quantile_for_mixed_species <- 1-quantiles_by_singlets$fraction_removed[quantiles_by_singlets$library=='63_20-hg19']
others <- ataqv[ataqv$library %in% MIXED_SPECIES_LIBRARIES & !ataqv$library %in% quantiles_by_singlets$library,]
thresholds_for_others <- others %>%
  dplyr::group_by(library) %>%
  dplyr::summarize(max_hqaa_threshold=quantile(hqaa, c(quantile_for_mixed_species)))
thresholds <- bind_rows(thresholds, thresholds_for_others)

# display as a stacked barplot...
atac$is_doublet <- ifelse(atac$assignment=='doublet', 'doublet', 'singlet')
p <- ggplot(atac) +
  geom_histogram(aes(x=hqaa, fill=is_doublet)) +
  facet_wrap(~library) +
  scale_x_log10() +
  scale_fill_viridis_d(option='C') +
  geom_vline(aes(xintercept=max_hqaa_threshold), data=thresholds[thresholds$library %in% atac$library,], linetype='dashed')
pdf('atac-doublet-stacked-barplots.pdf')
p
dev.off()


p <- ggplot(ataqv) +
  geom_histogram(aes(x=hqaa), bins=50) +
  facet_wrap(~library, scales='free_y') +
  scale_x_log10() +
  scale_fill_viridis_d(option='C') +
  geom_vline(aes(xintercept=max_hqaa_threshold), data=thresholds, linetype='dashed')
pdf('atac-max-hqaa-thresholds.pdf')
p
dev.off()


# write out the final list of nuclei
write.table(thresholds, file='atac-max-hqaa-thresholds.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = T)
atac_final_nuclei <- left_join(ataqv, thresholds) %>%
  dplyr::filter(hqaa<max_hqaa_threshold) %>%
  dplyr::select(library, barcode) %>%
  left_join(nuclei_assignments)
atac_final_nuclei$assignment[is.na(atac_final_nuclei$assignment) & grepl('hg19', atac_final_nuclei$library)] <- 'KSM1'
atac_final_nuclei$assignment[grepl('rn6', atac_final_nuclei$library)] <- 'rat1'
atac_final_nuclei <- atac_final_nuclei %>% dplyr::filter(assignment!='doublet')
atac_final_nuclei <- atac_final_nuclei[order(atac_final_nuclei$library, atac_final_nuclei$barcode),]
write.table(atac_final_nuclei, file='atac-nuclei-with-individuals.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = F)