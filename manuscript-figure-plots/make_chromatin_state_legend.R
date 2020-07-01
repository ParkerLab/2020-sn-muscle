#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)

args <- commandArgs(T)
chromhmm_gb_bed_file <- args[1]

chromhmm <- read.table(chromhmm_gb_bed_file, head = F, as.is = T, sep = '\t')
chromhmm <- unique(chromhmm[,c('V4', 'V9')])
colnames(chromhmm) <- c('state', 'rgb')
chromhmm$state <- gsub('_', ' ', chromhmm$state)
chromhmm$order <- as.numeric(gsub('(\\d+) .*', '\\1', chromhmm$state))
chromhmm$state <- gsub('\\d+ ', '', chromhmm$state)
chromhmm <- chromhmm[order(chromhmm$order),]
chromhmm$state <- factor(chromhmm$state, levels = chromhmm$state, ordered = T)
chromhmm$r <- as.numeric(unlist(lapply(strsplit(chromhmm$rgb, ','), function(x){x[1]})))
chromhmm$g <- as.numeric(unlist(lapply(strsplit(chromhmm$rgb, ','), function(x){x[2]})))
chromhmm$b <- as.numeric(unlist(lapply(strsplit(chromhmm$rgb, ','), function(x){x[3]})))
colors <- apply(chromhmm[,c('r', 'g', 'b')], 1 ,function(x){rgb(red = x[1], green = x[2], blue = x[3], maxColorValue = 255)})
names(colors) <- chromhmm$state
p <- ggplot(chromhmm) + geom_bar(aes(x = state, fill = state), color = 'black') +
  scale_fill_manual(values = colors) + guides(fill=guide_legend(title='Chromatin\nstate')) +
  theme(legend.position = 'top')
pdf('chrom_state_legend.pdf', height = 5, width = 7.5)
p
dev.off()
