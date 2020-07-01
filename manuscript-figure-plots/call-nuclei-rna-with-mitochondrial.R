#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(glue)
library(tidyr)
library(ggpointdensity)
library(viridis)

args <- commandArgs(T)
INITIAL_THRESHOLDS <- args[1]
LIBRARY_LABELS <- args[2]
QC <- args[3:length(args)]

qc <- bind_rows(lapply(QC, function(f){
  lib <- gsub('.txt', '', basename(f))
  tmp <- read.table(f, head = T, sep = '\t', as.is = T)
  tmp$library <- lib
  return(tmp)
}))

thresholds <- read.table(INITIAL_THRESHOLDS, head=T, sep='\t', stringsAsFactors = T)
library_labels <- read.table(LIBRARY_LABELS, head=T, sep='\t', stringsAsFactors = F)
library_labels <- library_labels[library_labels$modality=='RNA',]
qc <- qc[qc$number_umis >=10,]


# Plot for all libraries we used downstream...
used_downstream <- qc[qc$library %in% library_labels$library,]
used_downstream <- left_join(used_downstream, library_labels[,c('library', 'name')])
#used_downstream$name <- gsub('no FANS', 'crude', used_downstream$name)
used_downstream_thresholds <- thresholds[thresholds$library %in% library_labels$library,] %>%
  left_join(used_downstream %>% dplyr::select(name, library) %>% unique())

p <- ggplot(used_downstream) + geom_pointdensity(aes(x = number_umis, y=fraction_mitochondrial+0.001)) +
  geom_vline(aes(xintercept = min_hqaa), data = used_downstream_thresholds, color='red', linetype='dashed') +
  geom_vline(aes(xintercept = max_hqaa), data = used_downstream_thresholds, color='red', linetype='dashed') +
  geom_hline(aes(yintercept = max_mitochondrial+0.001), data = used_downstream_thresholds, color='red', linetype='dashed') +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~name) +
  theme_bw() +
  xlab('Number of UMIs') +
  ylab('Fraction mitochondrial + 0.001') +
  theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
  scale_color_viridis(trans='log') +
  guides(color=guide_colorbar(title='log(n_neighbors)'))
png('umis-vs-mitochondrial-used-downstream.png', height = 5, width = 8, res=300, units='in')
p
dev.off()



p <- ggplot(used_downstream) + 
  geom_histogram(aes(x = number_umis), bins=100) +
  geom_vline(aes(xintercept = min_hqaa), data = used_downstream_thresholds, color='red', linetype='dashed') +
  geom_vline(aes(xintercept = max_hqaa), data = used_downstream_thresholds, color='red', linetype='dashed') +
  scale_x_log10() +
  facet_wrap(~name, scales='free_y') +
  theme_bw() +
  xlab('Number of UMIs') +
  coord_cartesian(ylim=c(0, 5000))
png('umis-used-downstream.png', height = 5, width = 5, res=300, units='in')
p
dev.off()





# First, do the plots for the FANS vs no FANS
fans_vs_no_fans <- qc[qc$library %in% c('133155-hg19', '133156-hg19', '133157-hg19', '133158-hg19'),] %>%
  #dplyr::mutate(fans_status=ifelse(library %in% c('133155-hg19', '133157-hg19'), 'Crude', 'FANS'),
  dplyr::mutate(fans_status=ifelse(library %in% c('133155-hg19', '133157-hg19'), 'no FANS', 'FANS'),
                replicate=ifelse(library %in% c('133155-hg19', '133156-hg19'), 'rep. 1', 'rep. 2'))
fans_vs_no_fans_thresholds <- thresholds[thresholds$library %in% fans_vs_no_fans$library,] %>%
  left_join(fans_vs_no_fans %>% dplyr::select(fans_status, library, replicate) %>% unique())

p <- ggplot(fans_vs_no_fans) + geom_pointdensity(aes(x = number_umis, y=fraction_mitochondrial+0.001)) +
  geom_vline(aes(xintercept = min_hqaa), data = fans_vs_no_fans_thresholds, color='red', linetype='dashed') +
  geom_vline(aes(xintercept = max_hqaa), data = fans_vs_no_fans_thresholds, color='red', linetype='dashed') +
  geom_hline(aes(yintercept = max_mitochondrial+0.001), data = fans_vs_no_fans_thresholds, color='red', linetype='dashed') +
  scale_x_log10() +
  scale_y_log10() +
  facet_grid(replicate~fans_status) +
  theme_bw() +
  xlab('Number of UMIs') +
  ylab('Fraction mitochondrial + 0.001') +
  theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
  scale_color_viridis(trans='log') +
  guides(color=guide_colorbar(title='log(n_neighbors)'))
png('umis-vs-mitochondrial-fans-vs-no-fans.png', height = 5, width = 5, res=300, units='in')
p
dev.off()



p <- ggplot(fans_vs_no_fans) + 
  geom_histogram(aes(x = number_umis), bins=100) +
  geom_vline(aes(xintercept = min_hqaa), data = fans_vs_no_fans_thresholds, color='red', linetype='dashed') +
  geom_vline(aes(xintercept = max_hqaa), data = fans_vs_no_fans_thresholds, color='red', linetype='dashed') +
  scale_x_log10() +
  facet_grid(replicate~fans_status) +
  theme_bw() +
  xlab('Number of UMIs') +
  coord_cartesian(ylim=c(0, 5000))
png('umis-fans-vs-no-fans.png', height = 5, width = 5, res=300, units='in')
p
dev.off()


# Now, do the plots for 20k vs 40k
loading <- qc[qc$library %in% c('63_20_rna-hg19', '63_40_rna-hg19'),] %>%
  dplyr::mutate(nuclei_loaded=paste0(gsub('63_(\\d+)_rna-hg19', '\\1', library), 'k nuclei'))
loading_thresholds <- thresholds[thresholds$library %in% loading$library,] %>%
  left_join(loading %>% dplyr::select(nuclei_loaded, library) %>% unique())

p <- ggplot((loading)) + geom_pointdensity(aes(x = number_umis, y=fraction_mitochondrial+0.001)) +
  geom_vline(aes(xintercept = min_hqaa), data = loading_thresholds, color='red', linetype='dashed') +
  geom_vline(aes(xintercept = max_hqaa), data = loading_thresholds, color='red', linetype='dashed') +
  geom_hline(aes(yintercept = max_mitochondrial+0.001), data = loading_thresholds, color='red', linetype='dashed') +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~nuclei_loaded) +
  theme_bw() +
  xlab('Number of UMIs') +
  ylab('Fraction mitochondrial + 0.001') +
  theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
  scale_color_viridis(trans='log') +
  guides(color=guide_colorbar(title='log(n_neighbors)'))

png('umis-vs-mitochondrial-20k-vs-40k.png', height = 3, width = 5, res=300, units='in')
p
dev.off()


p <- ggplot(loading) + 
  geom_histogram(aes(x = number_umis), bins=100) +
  geom_vline(aes(xintercept = min_hqaa), data = loading_thresholds, color='red', linetype='dashed') +
  geom_vline(aes(xintercept = max_hqaa), data = loading_thresholds, color='red', linetype='dashed') +
  scale_x_log10() +
  facet_wrap(~nuclei_loaded, scales='free_y') +
  theme_bw() +
  xlab('Number of UMIs') +
  coord_cartesian(ylim=c(0, 5000))
png('umis-20k-vs-40k.png', height = 3, width = 5, res=300, units='in')
p
dev.off()



# also write out the table of thresholds used...
supplemental_table_of_thresholds_used <- left_join(library_labels, thresholds) %>%
  dplyr::select(library, experiment, min_hqaa, max_hqaa, max_mitochondrial)
supplemental_table_of_thresholds_used$max_mitochondrial <- supplemental_table_of_thresholds_used$max_mitochondrial - 0.001
colnames(supplemental_table_of_thresholds_used) <- c('library', 'experiment', 'min. UMIs', 'max. UMIs', 'max. fraction mitochondrial')
write.table(supplemental_table_of_thresholds_used, file = 'rna-qc-thresholds.csv', append=F, quote = F, sep = ',', row.names = F, col.names = T)
