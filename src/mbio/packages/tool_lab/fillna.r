# Title     : sub-script of table.kit.py
# Objective : use for filling NA values of data frame
# Created by: jincheng.qin
# Created on: 2019/7/24

# load package
library(getopt)
library(missForest)
library(impute)
library(pcaMethods)
library(pcaMethods)
# set options
command <- matrix(c(
  "input", "i", 1, "character", "text table file",
  "method", "m", 1, "character", "method for filling NA values",
  "output", "o", 1, "character", "NA filled table",
  "help", "h", 0, "logical", "show this help message and exit"
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

# check options
if (!is.null(opts$help)) {
  cat(getopt(command, usage = TRUE))
  q(status = -1)
}
if (is.null(opts$input) || is.null(opts$method) || is.null(opts$output)) {
  cat(getopt(command, usage = TRUE))
  q(status = -2)
}

# prepare data
print(paste("start reading", opts$input))
df.raw <- read.table(opts$input, header = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)
if (opts$method == 'missforest'){    # Deal with characters, added by zhangyitong on 20211202
    df.raw[,2:ncol(df.raw)] <- sapply(df.raw[,2:ncol(df.raw)], as.numeric, USE.NAMES=FALSE)
}
vec.cn <- colnames(df.raw)
vec.rn <- df.raw[,1]
df.int <- df.raw[,-1]

# main
print(paste("start filling NA values by", opts$method, "method"))
if (opts$method == 'missforest') {
  df.new <- missForest(df.int)$ximp
} else if (opts$method == "knn") {
  df.new <- t(impute::impute.knn(t(df.int))$data)
} else if (opts$method == "bpca") {
  df.new <- pca(df.int, nPcs = 5, method = "bpca", center = TRUE)@completeObs
} else if (opts$method == "ppca") {
  df.new <- pca(df.int, nPcs = 5, method = "ppca", center = TRUE)@completeObs
} else if (opts$method == "nipals") {
  df.new <- pca(df.int, nPcs = 5, method = "nipals", center = TRUE)@completeObs
}

# export result
print(paste("start exporting", opts$output))
df.out <- data.frame(vec.rn, df.new)
colnames(df.out) <- vec.cn
write.table(df.out, opts$output, quote = FALSE, sep = "\t", row.names = FALSE)
