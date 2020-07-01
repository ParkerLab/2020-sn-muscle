#!/usr/bin/env Rscript
library(ggplot2)
library(glue)
library(dplyr)
library(tidyr)

args <- commandArgs(T)

SOUPORCELL <- args[1]
DOUBLETFINDER <- args[2]

souporcell <-  read.table(SOUPORCELL, head = F, sep = '\t', col.names = c('library', 'barcode', 'assignment'), colClasses = c('character', 'character', 'character'))
doubletfinder <-  read.table(DOUBLETFINDER, head = F, sep = '\t', col.names = c('library', 'barcode', 'assignment'), colClasses = c('character', 'character', 'character'))

# compare demuxlet to souporcell results...
#both <- left_join(doubletfinder %>% dplyr::rename(df=assignment), souporcell %>% dplyr::rename(sc=assignment))
#both$sc[both$sc!='doublet' & !is.na(both$sc)] <- 'singlet'
#both <- both[!is.na(both$sc),]
# calculate concordance between DF and SC:
#concordance <- both %>% dplyr::group_by(library, df, sc) %>%
#  dplyr::summarize(count=n())
#print('Concordance:')
#print(concordance)

# any nucleus that has been assigned as a doublet by either demuxlet or souporcell
# should be considered a doublet...
both <- left_join(doubletfinder %>% dplyr::rename(df=assignment), souporcell %>% dplyr::rename(sc=assignment))
both$sc[is.na(both$sc)] <- 'KSM1' # these libraries were not mixed
both$assignment <- both$sc
both$assignment[both$df!='singlet'] <- 'doublet'

print('Final assignments:')
both %>%
  dplyr::group_by(library, df, sc, assignment) %>%
  dplyr::summarize(count=n())


counts <- both %>%
  dplyr::group_by(library, assignment) %>%
  dplyr::summarize(count=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(library) %>%
  dplyr::mutate(fraction=count/sum(count),
                perc=round(fraction*100, 1)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(label=glue('{count} ({perc}%)'))

#counts$library <- gsub('63_20_rna-hg19', '20k input', counts$library)
#counts$library <- gsub('63_40_rna-hg19', '40k input', counts$library)

p <- ggplot(counts) +
  geom_bar(aes(x=assignment, y=fraction, fill=library), position=position_dodge(width=1), stat='identity') +
  geom_text(aes(x=assignment, y=fraction+0.01, label=label, fill=library), position=position_dodge(width=1)) +
  theme_bw() +
  ylab('Fraction of pass-QC nuclei') +
  xlab('Sample or doublet') +
  guides(fill=guide_legend(title=''))
pdf('snRNA-doublets.pdf', height=5, width=3)
print(p)
dev.off()

write.table(both[,c('library', 'barcode', 'assignment')], file='rna-souporcell-assignments.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = F)
