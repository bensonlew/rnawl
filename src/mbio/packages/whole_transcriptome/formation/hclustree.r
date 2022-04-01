# accept parameters

library(getopt)

command <- matrix(c(
  'input', 'i', 1, 'character', 'input file containing expression matrix',
  'dist_method', 'd', 1, 'character', 'distance measure to be used',
  'clus_method', 'c', 1, 'character', 'agglomeration method to be used',
  'output', 'o', 1, 'character', 'output file containing tree string',
  'help', 'h', 0, 'logical', 'show this help message and exit'
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

if (!is.null(opts$help)) {
  cat(getopt(command, usage = TRUE))
  q(status=2)
}
if (is.null(opts$input) || is.null(opts$output) || is.null(opts$dist_method) || is.null(opts$clus_method)) {
  cat(getopt(command, usage = TRUE))
  q(status=2)
}

# load packages

library(ape)
library(fastcluster)

# main

print("start loading required data")
exp.df <- read.delim(opts$input, row.names = 1)
print("start computing distance matrix")
dat.dist <- dist(t(exp.df), method = opts$dist_method)
print("start implementing hierarchical clustering")
ret.hc <- hclust(dat.dist, method = opts$clus_method)
print("start converting hclust object to phylo object")
ret.phy <- as.phylo(ret.hc)
print("start writing a tree in parenthetic format")
write.tree(ret.phy, file = opts$output)
print(paste("succeed in exporting", opts$output))
