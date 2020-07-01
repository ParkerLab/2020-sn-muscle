#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)

args <- commandArgs(T)
intersected <- args[1]
roadmap_cell_types <- args[2]

intersected <- read.table(intersected, head=F, as.is=T, sep='\t', col.names = c('chrom', 'start', 'end', 'eid', 'state'))
roadmap_cell_types <- read.table(roadmap_cell_types, head=T, sep='\t', as.is=T, comment.char = '')
#roadmap_cell_types[roadmap_cell_types$comment!='',]
roadmap_cell_types <- roadmap_cell_types %>% dplyr::select(eid, group, anatomy, standardized_name, type)
intersected <- left_join(intersected, roadmap_cell_types)
print(intersected[intersected$anatomy=='MUSCLE',])
print(intersected[intersected$anatomy=='KIDNEY',])
print(intersected[grep('Enh', intersected$state),])
