#!/usr/bin/env Rscript

library(clusterProfiler)
library(enrichplot)
library(getopt)

command <- matrix(c(
  'geneList', 'g', 1, 'character', 'A csv file contains two columns, one for gene ID (no duplicated allowed) and another one for fold change.',
  'organism', 'o', 1, 'character', 'It can be any species that have KEGG annotation data available in KEGG database.',
  'nPerm', 'n', 2, 'integer', 'The number of permutations, default 1000.',
  'minGSSize', 'm', 2, 'integer', 'Gene sets smaller than this number are EXLCUDED from the analysis.',
  'pvalueCutoff', 'p', 2, 'double', 'pvalue cutoff on enrichment tests to report as significant.',
  'nGenesets', 's', 2, 'integer', 'The number of gene sets to be displayed on the same figure.',
  'help', 'h', 0, 'logical', 'show this help message and exit.'
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

if (!is.null(opts$help)) {
  cat(getopt(command, usage = TRUE))
  q(status=-1)
}
if (is.null(opts$geneList) || is.null(opts$organism)) {
  cat(getopt(command, usage = TRUE))
  q(status=-2)
}
## 设置默认值
if ( is.null(opts$nPerm) ) {
  opts$nPerm = 1000
}
if ( is.null(opts$pvalueCutoff) ) {
  opts$pvalueCutoff = 0.05
}
if ( is.null(opts$minGSSize) ) {
  opts$minGSSize = 10
}
if ( is.null(opts$nGenesets) ) {
  opts$nGenesets = 5
}

gene = read.table(opts$geneList)
geneList <- gene[,2]
names(geneList) <- as.character(gene[,1])
geneList <- sort(geneList, decreasing = TRUE)

result <- gseKEGG(geneList     = geneList,
               organism     = opts$organism,
               nPerm        = opts$nPerm,
               minGSSize    = opts$minGSSize,
               pvalueCutoff = opts$pvalueCutoff,
               verbose      = TRUE)

write.table(as.data.frame(result), "gsea.txt", row.names = F, sep="\t", quote=F)
pdf(file = 'gseaplot.pdf',height = 10,width = 10)
gseaplot2(result, geneSetID = 1:opts$nGenesets, pvalue_table=TRUE)
dev.off()
pdf(file = 'gseaupset.pdf',height = 10,width = 10)
upsetplot(result, n = opts$nGenesets)
dev.off()