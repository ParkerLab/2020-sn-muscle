#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
make_option(c("--rda"), action = "store", type = "character", help = "[Required] rda file from liger_factorize"),
make_option(c("--prefix"), action = "store", type = "character", default = '', help = "[Required] Prefix of the output files.")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts <- parse_args(option_parser)

# for testing
# opts <- list(
#   'rda' = '/lab/work/porchard/sn-muscle-project/work/downstream/results/liger/round-1/reclustered.Rda',
#   'prefix' = 'test'
# )

PREFIX <- opts$prefix

library(liger)
library(glue)

load(opts$rda)

pdf(glue('{PREFIX}plot-gene-loadings.pdf'))
plotGeneLoadings(int.muscle, return.plots = F)
dev.off()

pdf(glue('{PREFIX}plot-factors.pdf'))
plotFactors(int.muscle, num.genes = 10)
dev.off()

pdf(glue('{PREFIX}plot-cluster-proportions.pdf'))
plotClusterProportions(int.muscle)
dev.off()

pdf(glue('{PREFIX}plot-cluster-factors.pdf'))
plotClusterFactors(int.muscle, use.aligned=T)
dev.off()
