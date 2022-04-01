#!/usr/bin/env Rscript

library(edgeR)
library(NOISeq)
library(limma)
library(DESeq2)
library(getopt)

command <- matrix(c(
  'input', 'i', 1, 'character', 'input file exp_matrix',
  'method', 'm', 1, 'character', 'normalized method to be used',
  'output', 'o', 1, 'character', 'output file exp_normalized_matrix',
  'number', 'n', 1, 'numeric', 'Retention significant number',
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

# prepare data
print(paste("start reading", opts$input))
df.raw <- read.table(opts$input, header = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)
vec.cn <- colnames(df.raw)
vec.rn <- df.raw[,1]
df.int <- df.raw[,-1]
print(vec.cn)
# main
print(paste("start filling NA values by", opts$method, "method"))
#TMM
if (opts$method == 'TMM') {
exp_factor <- calcNormFactors(df.int, method = 'TMM')
exp_norm <- t(t(df.int)/exp_factor)
}
#TMMwzp
if (opts$method == 'TMMwzp') {
exp_factor <- calcNormFactors(df.int, method = 'TMMwzp')
exp_norm <- t(t(df.int)/exp_factor)
}
#Upper Quartile
if (opts$method == 'uqua') {
exp_factor <- calcNormFactors(df.int, method = 'upperquartile')
exp_norm <- t(t(df.int)/exp_factor)
}
#RLE
if (opts$method == 'RLE') {
exp_factor <- calcNormFactors(df.int, method = 'upperquartile')
exp_norm <- t(t(df.int)/exp_factor)
}
#DESeq2
if (opts$method == 'DESeq2') {
print(ncol(df.int))

condition <- factor(c(rep(1,2),rep(2,length(vec.cn)-1-2)))
print(condition)
df.int <- floor(df.int)
dds <- DESeqDataSetFromMatrix(df.int, DataFrame(condition), ~ condition)
dds <- DESeq(dds)
exp_norm <- counts(dds, normalized=TRUE)

}



# export result
print(paste("start exporting", opts$output))
exp_norm <- round(exp_norm, opts$number)
df.out <- data.frame(vec.rn, exp_norm)
colnames(df.out) <- vec.cn
write.table(df.out, opts$output, quote = FALSE, sep = "\t", row.names = FALSE)


