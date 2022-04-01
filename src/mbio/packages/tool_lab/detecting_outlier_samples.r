# 安装步骤如下：
# export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/rna/miniconda3/bin:$PATH
# install.packages("flock")
# install.packages("RhpcBLASctl")
# install.packages("bigparallelr")
# install.packages("nabor")
# install.packages("robustbase")
# install.packages("RSpectra")
# install.packages("bigstatsr")
# install.packages("/mnt/ilustre/users/sanger-dev/app/install_packages/privefl-bigutilsr-069e8f0.tar.gz",repo=NULL,type="source")

library(bigutilsr)
library(ggplot2)
library(bigstatsr)
library(magrittr)
exp <- read.table("/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/tmp/ref.gene.tpm.matrix", header=T, sep="\t", row.names=1)
exp_matrix <- data.matrix(exp)
exp_matrix <- exp_matrix[which(rowSums(exp_matrix) > 0),]
exp <- t(exp_matrix)
pca <- prcomp(exp, scale. = TRUE, rank. = 10)
U <- pca$x
theme_set(bigstatsr::theme_bigstatsr(0.8))

ind.out <- apply(U, 2, function(x) which( (abs(x - median(x)) / mad(x)) > 6 )) %>%
  Reduce(union, .) %>%
  print()
col <- rep("black", nrow(U)); col[ind.out] <- "red"
d = qplot(U[, 1], U[, 2], color = I(col), size = I(1)) + coord_equal() + geom_text(aes(label = rownames(U)),cex=1)
ggsave("outlier.pdf", d, width = 6, height = 6)
dist <- apply(U, 2, function(x) abs(x - median(x)) / mad(x)) %>%
  apply(1, max)
qplot(U[, 1], U[, 3], color = dist, size = I(3)) + coord_equal() +
  scale_color_viridis_c(trans = "log", breaks = c(1, 3, 6))
qplot(y = sort(dist, decreasing = TRUE)) +
  geom_hline(yintercept = 6, color = "red")
