#!/usr/bin/env Rscript
library(ggplot2)
library(getopt)
library(ggrepel)

command <- matrix(c(
  'diff_file', 'g', 1, 'character', 'A csv file contains three columns, column one for gene ID (no duplicated allowed), column two for fold change and column for pvalue.',
  "fc", "f", 2, "numeric", "log2fc",
  "top", "t", 2, "integer", "top gene to add label",
  "method", "m", 1, "integer", "color schemes",
  "x_axis_name", "x", 2, "character", "x_axis_name",
  "y_axis_name", "y", 2, "character", "y_axis_name",
  'help', 'h', 0, 'logical', 'show this help message and exit.'
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

if (!is.null(opts$help)) {
  cat(getopt(command, usage = TRUE))
  q(status=-1)
}
if (is.null(opts$diff_file)) {
  cat(getopt(command, usage = TRUE))
  q(status=-2)
}
## 设置默认值
if ( is.null(opts$x_axis_name) ) {
  opts$x_axis_name = "Rank of differentially expressed genes"
}
if ( is.null(opts$y_axis_name) ) {
  opts$y_axis_name = "Log2FoldChange"
}
if ( is.null(opts$fc) ) {
  opts$fc = 1
}
if ( is.null(opts$top) ) {
  opts$top = 5
}
if ( is.null(opts$method) ) {
  opts$method = 1
}
if (opts$method == 1 ){
    colo_low = "red"
    colo_mid = "grey"
    colo_up = "blue"
   }else if (opts$method == 2){
    colo_low = "red"
    colo_mid = "grey"
    colo_up = "green"
   }else if (opts$method == 3){
    colo_low = "red"
    colo_mid = "black"
    colo_up = "blue"
   }else if(opts$method == 4){
    colo_low = "red"
    colo_mid = "black"
    colo_up = "green"
   }

diff_file = read.table(opts$diff_file, header=T)
colnames(diff_file) <- c("gene_id", "log2fc", "pvalue")
diff_file_sort = diff_file[order(diff_file$log2fc,decreasing = TRUE),]
diff_file_sort$rank = c(1:length(diff_file_sort$gene_id))
diff_file_sort$absfc = abs(diff_file_sort$log2fc)
label_df = merge(head(diff_file_sort, n=opts$top), tail(diff_file_sort, n=opts$top), all=TRUE)
pdf(file = 'diff_plot.pdf',height = 5,width = 5)
ggplot(diff_file_sort, aes(x=rank, y=log2fc)) +
    geom_point(aes(color = pvalue, size=absfc)) +
    scale_size_continuous(range = c(0.1,3)) +
    scale_color_gradient2(low = colo_low, mid = colo_mid, high = colo_up, midpoint = 0.5) +
    labs(x = opts$x_axis_name, y = opts$y_axis_name, color = 'Pvalue', size = '|Log2FC|') +
    geom_hline(aes(yintercept = opts$fc), linetype = 2) +
    geom_hline(aes(yintercept = -opts$fc), linetype = 2) +
    geom_hline(yintercept = 0) +
    geom_vline(aes(xintercept = floor(length(gene_id)/2)), linetype = 2) +
    geom_text_repel(label_df, mapping = aes(rank, log2fc, label = gene_id, size = 8), max.overlaps=opts$top*2)
dev.off()