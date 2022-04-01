#!/usr/bin/env Rscript

library(getopt)

command <- matrix(c(
  'pvalue', 'p', 1, 'character', 'numeric vector of p-values (possibly with NAs). Any other R is coerced by as.numeric.',
  'method', 'm', 1, 'character', 'correction method. Can be abbreviated.',
  'help', 'h', 0, 'logical', 'show this help message and exit'
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

if (!is.null(opts$help)) {
  cat(getopt(command, usage = TRUE))
  q(status=-1)
}

pvalue = read.table(opts$p)
p = as.numeric(pvalue[, 1])
padjust = p.adjust(p, method = opts$m, n = length(p))
result = data.frame(pvalue=p, padjust=padjust)
col_names = c("pvalue", opts$m)
colnames(result) <- col_names
output = paste(opts$m, "txt", sep=".")
write.table(result, output, row.names = F, sep="\t", quote=F)
