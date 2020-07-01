#!/usr/bin/env Rscript
#library(optparse)

#option_list <- list(
#  make_option(c("--out"), action = "store", type = "character", help = "[Required] Name of output pdf"),
#  make_option(c("--height"), action = "store", type = "numeric", default = 15, help = "[Optional] Height (in inches) of output pdf"),
#  make_option(c("--width"), action = "store", type = "numeric", default = 15, help = "[Optional] Width (in inches) of output pdf")
#)

#option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
#opts <- parse_args(option_parser)

library(dplyr)
library(ggplot2)
library(glue)

args <- commandArgs(T)
LIBRARY_LABELS <- args[1]
OVERLAP_FILES <- args[2:length(args)]

library_labels <- read.table(LIBRARY_LABELS, head = T, as.is = T, sep = '\t')
#library_labels$fans_status <- ifelse(grepl('no FANS', library_labels$name), 'Crude', 'FANS')
library_labels$fans_status <- ifelse(grepl('no FANS', library_labels$name), 'no FANS', 'FANS')
library_labels$loading_concentration <- NA
library_labels$loading_concentration[grepl('20', library_labels$name)] <- '20k'
library_labels$loading_concentration[grepl('40', library_labels$name)] <- '40k'
  
  
overlap <- bind_rows(lapply(OVERLAP_FILES, function(f){
  tmp <- read.table(f, head = T, as.is = T, sep = '\t')
  tmp$library <- gsub('.chromhmm_overlap.txt', '', basename(f))
  return(tmp)
}))

overlap <- overlap %>%
  dplyr::group_by(library, tissue, tss_relative) %>%
  dplyr::mutate(overlap=overlap/sum(overlap)) %>%
  dplyr::ungroup()


unique_types <- unique(overlap$state)
factor_order <- rev(c("Active TSS", "Weak TSS", "Flanking TSS", "Bivalent poised TSS", "Genic enhancer", "Active enhancer 1", "Active enhancer 2", "Weak enhancer", "Strong transcription", "Weak transcription", "Repressed polycomb", "Weak repressed polycomb", "Quiescent low signal"))
overlap$state <- as.vector(sapply(overlap$state, function(x){paste(strsplit(x, "_")[[1]][-1], collapse=" ")}))
overlap$state <- factor(overlap$state, levels=factor_order)
overlap <- overlap[order(overlap$state),]

color_pallette <- rev(c("#F91000", "#FA5D5F", "#FA5D5F", "#C90FBE", "#F9C202", "#F9C202", "#F9C202", "#FFFB04", "#24A647", "#8FFF5C", "#757575", "#B8B8B8", "#FFFFFF"))

overlap$tss_relative <- paste('TSS', overlap$tss_relative, 'peaks')

# plot for FANS...
fans <- left_join(overlap, library_labels) %>%
  dplyr::filter(!is.na(experiment) & experiment=='FANS vs. no FANS') %>%
  dplyr::mutate(label=glue('{fans_status}, {replicate}')) %>%
  dplyr::select(-experiment, -modality, -species) %>%
  dplyr::filter(tss_relative=='TSS distal peaks')

p <- ggplot(fans) + 
  geom_bar(aes(x=tissue,y=overlap, fill=state), stat="identity", colour="black") + 
  coord_flip() + 
  facet_grid(replicate ~ fans_status) + 
  scale_fill_manual(values=color_pallette) +  
  theme_grey(base_size = 23) + 
  ylab("Proportion overlap with chromatin state") + 
  xlab("Tissue") +
  guides(fill=guide_legend(title='Chromatin state'))
pdf('fans-chromhmm-overlap.pdf', width=15, height=13)
print(p)
dev.off()

# plot for 20k vs 40k
loading <- left_join(overlap, library_labels) %>%
  dplyr::filter(!is.na(experiment) & experiment=='20k vs. 40k') %>%
  dplyr::mutate(label=glue('{loading_concentration}')) %>%
  dplyr::select(-experiment, -modality, -species) %>%
  dplyr::filter(tss_relative=='TSS distal peaks')

p <- ggplot(loading) + 
  geom_bar(aes(x=tissue,y=overlap, fill=state), stat="identity", colour="black") + 
  coord_flip() + 
  facet_grid(. ~ label) + 
  scale_fill_manual(values=color_pallette) +  
  theme_grey(base_size = 23) + 
  ylab("Proportion overlap with chromatin state") + 
  xlab("Tissue") +
  guides(fill=guide_legend(title='Chromatin state'))
pdf('loading-chromhmm-overlap.pdf', width=15, height=10)
print(p)
dev.off()
