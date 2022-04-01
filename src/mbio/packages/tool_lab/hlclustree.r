#!/usr/bin/env Rscript

library(ape)
library(getopt)

command <- matrix(c(
  'input', 'i', 1, 'character', 'input file containing distance data between samples',
  'method', 'm', 1, 'character', 'agglomeration method to be used',
  'output', 'o', 1, 'character', 'output file containing tree string',
  'help', 'h', 0, 'logical', 'show this help message and exit'
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

if (!is.null(opts$help)) {
  cat(getopt(command, usage = TRUE))
  q(status=-1)
}
if (is.null(opts$input) || is.null(opts$output) || is.null(opts$method)) {
  cat(getopt(command, usage = TRUE))
  q(status=-2)
}

dist.df <- read.delim(opts$input, row.names = 1)
dat.dist <- as.dist(dist.df)
ret.hc <- hclust(dat.dist, method = opts$method)
ret.phy <- as.phylo(ret.hc)
write.tree(ret.phy, file = opts$output)
