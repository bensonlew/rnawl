#!/usr/bin/env Rscript

library(clusterProfiler)
library(getopt)

command <- matrix(c(
  'gene', 'g', 1, 'character', 'a vector of gene id',
  'pvalueCutoff', 'p', 2, 'double', 'pvalue cutoff on enrichment tests to report',
  'pAdjustMethod', 'd', 2, 'character', 'one of holm, hochberg, hommel, bonferroni, BH, BY, fdr, none',
  'universe', 'u', 2, 'character', 'background genes. If missing, the all genes listed in the database (eg TERM2GENE
table) will be used as background.',
  'minGSSize', 'm', 2, 'integer', 'minimal size of genes annotated for testing',
  'maxGSSize', 'n', 2, 'integer', 'maximal size of genes annotated for testing',
  'qvalueCutoff', 'q', 2, 'double', 'qvalue cutoff on enrichment tests to report as significant. Tests must pass i)
pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on adjusted pvalues
and iii) qvalueCutoff on qvalues to be reported.',
  'TERM2GENE', 't', 1, 'character', 'user input annotation of TERM TO GENE mapping, a data.frame of 2 column
with term and gene',
  'TERM2NAME', 'r', 2, 'character', 'user input of TERM TO NAME mapping, a data.frame of 2 column with term
and name',
  'help', 'h', 0, 'logical', 'show this help message and exit'
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

if (!is.null(opts$help)) {
  cat(getopt(command, usage = TRUE))
  q(status=-1)
}
if (is.null(opts$gene) || is.null(opts$TERM2GENE)) {
  cat(getopt(command, usage = TRUE))
  q(status=-2)
}
## 设置默认值
if ( is.null(opts$pvalueCutoff) ) {
  opts$pvalueCutoff = 1
}
if ( is.null(opts$pAdjustMethod) ) {
  opts$pAdjustMethod = "BH"
}
if ( is.null(opts$minGSSize) ) {
  opts$minGSSize = 10
}
if ( is.null(opts$maxGSSize) ) {
  opts$maxGSSize = 500
}
if ( is.null(opts$qvalueCutoff) ) {
  opts$qvalueCutoff = 1
}

gene = read.table(opts$gene)
TERM2GENE = read.table(opts$TERM2GENE, sep="\t")
if ( ! is.null(opts$TERM2NAME) ) {
  TERM2NAME = read.table(opts$TERM2NAME, sep="\t")
}else{
  TERM2NAME = NULL
}

if ( ! is.null(opts$universe) ) {
  universe = as.character(read.table(opts$universe)[, 1])
}else{
  universe = NULL
}

result = enricher(gene[, 1], TERM2GENE=TERM2GENE, TERM2NAME=TERM2NAME, universe=universe, pvalueCutoff=opts$pvalueCutoff,
                  pAdjustMethod=opts$pAdjustMethod, minGSSize=opts$minGSSize, maxGSSize=opts$maxGSSize,
                  qvalueCutoff=opts$qvalueCutoff)
write.table(as.data.frame(result), "enrichment.txt", row.names = F, sep="\t", quote=F)
