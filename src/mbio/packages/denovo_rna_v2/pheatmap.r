# Title     : sub-script of table.kit.py
# Objective : use for filling NA values of data frame
# Created by: jincheng.qin
# Created on: 2019/7/24

# load package
library(getopt)
library(pheatmap)
library(ape)
# set options
command <- matrix(c(
  "exp_matrix", "i", 1, "character", "text exp_matrix with scale",
  "sample_tree", "s", 1, "character", "sample_tree",
  "seq_tree", "m", 1, "character", "seq_tree",
  "output", "o", 1, "character", "output",
  'group2sample', 'g', 1, "character", 'group2sample',
  'group2seq', 'q', 1, 'character', 'group2seq',
  "help", "h", 0, "logical", "show this help message and exit"
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

# check options
if (!is.null(opts$help)) {
  cat(getopt(command, usage = TRUE))
  q(status = -1)
}

#if (is.null(opts$input) || is.null(opts$method) || is.null(opts$output)) {
#  cat(getopt(command, usage = TRUE))
#  q(status = -2)
#}

exp <- read.table(opts$exp_matrix, header = T, row.names = 1, sep='\t')
if (opts$sample_tree == 'False'){
hc_sample <-FALSE
group_sample <-NA
}else{
sample_tree <- read.tree(opts$sample_tree)
hc_sample <- as.hclust(sample_tree)
group_sample <-read.table(opts$group2sample, row.names=1)
sample_tree_index <- sample_tree$tip.label
exp <- exp[, sample_tree_index]
}
if (opts$seq_tree == 'False'){
hc_seq <-FALSE
group_seq <-NA
} else if (opts$seq_tree == 'kmeans'){

hc_seq <- FALSE
group_seq <- read.table(opts$group2seq, row.name=1, header = T)
seq_tree_index <- row.names(group_seq)
exp <- exp[seq_tree_index, ]
group_seq <- NA
} else{
seq_tree <- read.tree(opts$seq_tree)
hc_seq <- as.hclust(seq_tree)
group_seq <- read.table(opts$group2seq, row.name=NULL, header = T)
seq_id <- group_seq$seq_id
group_seq <- data.frame(group2 = factor(group_seq$group2))
rownames(group_seq) = seq_id
seq_tree_index <- seq_tree$tip.label
exp <- exp[seq_tree_index, ]
}

#sample_tree <- read.tree(opts$sample_tree)

#seq_tree <- read.tree(opts$seq_tree)
#group_sample <-read.table(opts$group2sample, row.names=1)
#group_seq <- read.table(opts$group2seq, row.name=1, header = T)
#hc_sample <- as.hclust(sample_tree)
#hc_seq <- as.hclust(seq_tree)
#sample_tree_index <- sample_tree$tip.label
#seq_tree_index <- seq_tree$tip.label



#exp_tree <- exp[seq_tree_index, sample_tree_index]
pdf=paste(opts$output,"/heatmap",".pdf",sep="")
pdf(pdf,width=8,height=6)
par(mar=c(3,2,2,5))
pheatmap(exp, scale = 'none', cluster_cols=hc_sample, cluster_rows=hc_seq, color=colorRampPalette(colors=c("blue","white","red"))(100),show_rownames=F, annotation_col=group_sample, annotation_row=group_seq, annotation_legend = FALSE, annotation_names_col = FALSE, annotation_names_row = FALSE)
dev.off()


