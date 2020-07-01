#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option(c("--mats"), action = "store", type = "character", default = "", help = "[Required] (comma-separated) paths to the *.txt files for input (columns are nuclei, rows are genes"),
  make_option(c("--factorization_k"), action = "store", type = "numeric", default = 20, help = "[Required] k parameter for factorization"),
  make_option(c("--factorization_lambda"), action = "store", type = "numeric", default = 5, help = "[Required] lambda parameter for factorization"),
  make_option(c("--out"), action = "store", type = "character", default = '', help = "[Required] Output .Rda file.")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts <- parse_args(option_parser)

# for testing
# opts <- list(
#   'mats' = paste(list.files('/lab/work/porchard/sn-muscle-project/work/downstream-final/results/liger/round-1/input-merged/', full.names=T), collapse=','),
#   'factorization_k' = 25,
#   'factorization_lambda' = 5,
#   'prefix' = 'test'
# )

PREFIX <- opts$prefix

count_files <- sort(unlist(strsplit(opts$mats, ',')))

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(liger)
library(glue)


liger_in <- lapply(count_files, function(f){
	x <- read.table(f, sep='\t', as.is = T, head = T)
	rownames(x) <- x[,1]
	x <- x[,2:ncol(x)]
	x <- as.matrix(x)
	colnames(x) <- gsub('^X', '', colnames(x))
	colnames(x) <- gsub('\\.', '-', colnames(x))
	return(x)
})

names(liger_in) <- sapply(count_files, function(f){gsub('.txt', '', basename(f))})

print('Selecting genes with:')
SELECT_GENES <- 'RNA-KSM2.liger'
print(SELECT_GENES)

int.muscle <- createLiger(liger_in)
int.muscle <- normalize(int.muscle)
int.muscle <- selectGenes(int.muscle, datasets.use = grep(SELECT_GENES, names(liger_in))) # use RNA to select genes and ignore ATAC
int.muscle <- scaleNotCenter(int.muscle)

# factorization
print('Factorizing')
int.muscle <- optimizeALS(int.muscle, k = opts$factorization_k, lambda = opts$factorization_lambda, nrep=5)
print('Finished factorizing')
save(int.muscle, file=opts$out)
print('Finished saving object')
