#!/usr/bin/env Rscript
library(tidyr)
library(dplyr)
library(ggplot2)
library(glue)

args <- commandArgs(T)
CLUSTER_NAMES <- args[1]
ALL_MODELS <- args[2:length(args)]

cluster_names <- read.table(CLUSTER_NAMES, head=T, as.is=T, sep='\t')

load_result_file <- function(f) {
  tmp <- read.table(f, head = T, as.is = T, sep = '\t')
  NAME_RE <- '^(.*?)\\.(.*)\\.results$'
  tmp$model <- gsub(NAME_RE, '\\1', basename(f))
  tmp$trait <- gsub(NAME_RE, '\\2', basename(f))
  return(tmp)
}

results <- bind_rows(lapply(ALL_MODELS, load_result_file))
results <- results[grep('L2_1', results$Category),]
results$Category <- gsub('L2_1', '', results$Category)
results <- results[results$model=='tissues' | results$Category==results$model,]

for(i in 1:nrow(cluster_names)) {
  cl <- cluster_names$old_name[i]
  results$Category[results$Category==glue('cluster_{cl}')] <- cluster_names$new_name[i]
}

results$model[results$model!='tissues'] <- 'Joint model with open chromatin'
results$model[results$model=='tissues'] <- 'Joint model with open chromatin and all other cell types'
results$Category[results$Category=='liver'] <- 'Liver'
results$Category[results$Category=='adipose'] <- 'Adipose'
results$Category[results$Category=='beta_cell'] <- 'Beta cell'
results$Category <- factor(results$Category, levels=rev(c(cluster_names$new_name, unique(results$Category)[!unique(results$Category) %in% cluster_names$new_name])), ordered=T)

COLORS <- c('Joint model with open chromatin and all other cell types'='blue', 'Joint model with open chromatin'='red')

p <- ggplot(results) +
  geom_errorbar(aes(x = Category, color=model, ymin=Coefficient-1.96*Coefficient_std_error, ymax=Coefficient+1.96*Coefficient_std_error), position='dodge') +
  geom_point(aes(x = Category, color=model, y=Coefficient), position=position_dodge(0.9)) +
  theme_bw() +
  xlab('Cell type') + ylab('LDSC Coefficient') +
  coord_flip() +
  geom_hline(yintercept = 0, linetype='dashed') +
  scale_color_manual(values=COLORS) +
  facet_wrap(~trait, scales='free_x')
pdf(glue('all_traits.pdf'), width=15, height=15)
p
dev.off()


newnames <- c('diamante.T2Dbmiadj.European.ldsc'='T2D', 'manning_finsBMIadj'='Fasting insulin')
results <- results[results$trait %in% names(newnames),]
results$trait <- sapply(results$trait, function(x){newnames[x]})
results <- results[results$Category!='common_open_chromatin',]
# calculate significance
#number_tests <- nrow(results)
# adjust for
# 2 traits, (7 + 3) cell types across rat+human, 2 models (joint and single cell type)
number_tests <- 2 * 10 * 2
results$pvalue <- sapply(results$Coefficient_z.score, function(x){min(c(1, 2*pnorm(x, lower.tail = T), 2*pnorm(x, lower.tail = F)))})
results$bonferroni_significant <- results$pvalue <= (0.05/number_tests)

p <- ggplot(results) +
  geom_errorbar(aes(x = Category, color=model, ymin=Coefficient-1.96*Coefficient_std_error, ymax=Coefficient+1.96*Coefficient_std_error), position='dodge') +
  geom_point(aes(x = Category, color=model, y=Coefficient), position=position_dodge(0.9)) +
  geom_text(aes(x = Category, label='*', color=model, y=Coefficient+3*Coefficient_std_error), data=results[results$bonferroni_significant,], position=position_dodge(0.9)) +
  theme_bw() +
  xlab('') + ylab('LDSC Coefficient') +
  coord_flip() +
  geom_hline(yintercept = 0, linetype='dashed') +
  scale_color_manual(values=COLORS) +
  facet_wrap(~trait, scales='free_x') +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position = 'top', legend.direction = 'vertical') +
  guides(color=guide_legend(title=''))
pdf(glue('T2D-FIns.pdf'), width=6, height=4)
p
dev.off()



