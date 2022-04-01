#!/usr/bin/env Rscript

# load package
library(getopt)
# library(pheatmap)
library(circlize)
library(ComplexHeatmap)

# set options
command <- matrix(c(
  "data", "i", 1, "character", "input file containing data for drawing heatmap",
  "design", "d", 1, "character", "input file describing experimental design",
  "group", "n", 1, "character", "which group does the input data belong",
  "type", "t", 1, "character", "expression type of input data",
  "output", "o", 1, "character", "output directory for results",
  "help", "h", 0, "logical", "show this help message and exit"
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

# check options
if (!is.null(opts$help)) {
  cat(getopt(command, usage = TRUE))
  q(status=-1)
}
if (is.null(opts$data) || is.null(opts$design)) {
  cat(getopt(command, usage = TRUE))
  q(status=-2)
}

# main
print("start reading data table")
df.data <- read.delim(opts$data, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
print("start reading design table")
df.design <- read.delim(opts$design, header = TRUE, stringsAsFactors = FALSE)

print("start preparing arguments for plotting")
col.Time <- colorRamp2(
  c(min(df.design[["Time"]]), median(df.design[["Time"]]), max(df.design[["Time"]])),
  topo.colors(3)
)

int.unit <- max(apply(df.design['Sample'], 1, nchar)) + 6
top_annotation <- HeatmapAnnotation(
  Name = anno_text(df.design[["Sample"]], gp = gpar(fontface="bold"), rot = 90, offset = unit(int.unit, "mm")),
  Time = df.design[["Time"]],
  col = list(Time = col.Time),
  gp = gpar(col = "black"),
  show_legend = FALSE
)

heatmap_legend_param <- list(
  title = paste(opts$type, "(max scaled)", sep = "\n"),
  at = c(0, 0.25, 0.5, 0.75, 1),
  labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
  title_position = "topcenter"
)

col.Heatmap <- colorRamp2(c(0, 0.5, 1), c("blue", "yellow", "red"))

ht <- Heatmap(df.data, col = col.Heatmap,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = top_annotation,
  heatmap_legend_param = heatmap_legend_param
)

print("start exporting figures")
svg(paste(opts$output, paste("heatmap", opts$group, "svg", sep = "."), sep = "/"), width = 7, height = 14)
draw(ht, padding = unit(c(2, 2, int.unit + 6, 2), 'mm'))
# pheatmap(df.data, cluster_rows = FALSE, cluster_cols = FALSE,
         # annotation_col = df.design["Time"],
         # annotation_colors = list(Time = c("aquamarine", "brown")),
         # annotation_legend = FALSE,
         # show_rownames = FALSE)
dev.off()
png(paste(opts$output, paste("heatmap", opts$group, "png", sep = "."), sep = "/"), width = 4200, height = 8400, res = 600)
draw(ht, padding = unit(c(2, 2, int.unit + 6, 2), 'mm'))
# pheatmap(df.data, cluster_rows = FALSE, cluster_cols = FALSE,
         # annotation_col = df.design["Time"],
         # annotation_colors = list(Time = c("aquamarine", "brown")),
         # annotation_legend = FALSE,
         # show_rownames = FALSE)
dev.off()
pdf(paste(opts$output, paste("heatmap", opts$group, "pdf", sep = "."), sep = "/"), width = 7, height = 14, onefile = FALSE)
draw(ht, padding = unit(c(2, 2, int.unit + 6, 2), 'mm'))
# pheatmap(df.data, cluster_rows = FALSE, cluster_cols = FALSE,
         # annotation_col = df.design["Time"],
         # annotation_colors = list(Time = c("aquamarine", "brown")),
         # annotation_legend = FALSE,
         # show_rownames = FALSE)
dev.off()
