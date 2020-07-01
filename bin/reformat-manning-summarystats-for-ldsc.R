#!/usr/bin/env 
library(dplyr)
library(tidyr)

# reformat Manning 2012 SNPs for LDSC
fglu <- read.table(gzfile('MAGIC_Manning_et_al_FastingGlucose_MainEffect.txt.gz'), head=T, as.is=T, sep = '\t')
fglu_bmiadj <- fglu[,c('Snp', 'effect_allele', 'other_allele', 'maf', 'BMIadjMainEffects', 'BMIadjMainSE', 'BMIadjMainP')]
fglu_bmiadj <- fglu_bmiadj %>%
  dplyr::rename(effect=BMIadjMainEffects,
                SE=BMIadjMainSE,
                pvalue=BMIadjMainP)
fglu <- fglu[,c('Snp', 'effect_allele', 'other_allele', 'maf', 'MainEffects', 'MainSE', 'MainP')]
fglu <- fglu %>%
  dplyr::rename(effect=MainEffects,
                SE=MainSE,
                pvalue=MainP)
fglu$N <- 58074
fglu_bmiadj$N <- 58074
write.table(fglu_bmiadj, file='manning_fgluBMIadj.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = T)
write.table(fglu, file='manning_fglu.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = T)


fins <- read.table(gzfile('MAGIC_Manning_et_al_lnFastingInsulin_MainEffect.txt.gz'), head=T, as.is=T, sep = '\t')
fins_bmiadj <- fins[,c('Snp', 'effect_allele', 'other_allele', 'maf', 'BMIadjMainEffects', 'BMIadjMainSE', 'BMIadjMainP')]
fins_bmiadj <- fins_bmiadj %>%
  dplyr::rename(effect=BMIadjMainEffects,
                SE=BMIadjMainSE,
                pvalue=BMIadjMainP)
fins <- fins[,c('Snp', 'effect_allele', 'other_allele', 'maf', 'MainEffects', 'MainSE', 'MainP')]
fins <- fins %>%
  dplyr::rename(effect=MainEffects,
                SE=MainSE,
                pvalue=MainP)
fins$N <- 51750
fins_bmiadj$N <- 51750

write.table(fins_bmiadj, file='manning_finsBMIadj.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = T)
write.table(fins, file='manning_fins.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = T)
