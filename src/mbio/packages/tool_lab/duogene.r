#!/usr/bin/env Rscript

# load package
library(getopt)
library(ggplot2)
library(gggenes)

# set options
command <- matrix(c(
  "data", "i", 1, "character", "input file containing processed data",
  "output", "o", 1, "character", "output file path",
  "help", "h", 0, "logical", "show this help message and exit"
), byrow = TRUE, ncol = 5)
opts <- getopt(command)
example_genes = read.table(opts$data, sep='\t', header=TRUE)
pdf(file = opts$output,height = 8,width = 8)
ggplot(example_genes, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()
dev.off()



