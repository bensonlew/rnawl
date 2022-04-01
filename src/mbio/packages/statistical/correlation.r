library(pheatmap)
library(psych)
library(ape)

method="${method}" # pearson kendall
data=read.table("${inputfile}",header=T,row.names=1,sep="\t")

bin <- cor(data,method=method)
cor <- as.data.frame(bin)
#cor <- as.data.frame(bin$r)  #commented out code by khl 20170622
#t <- as.data.frame(bin$t)
#p <- as.data.frame(bin$p)
#cor <- round(cor,3)
write.table(cor,"${corr_matrix}",sep="\t", quote = F)
#write.table(p,"${pvalue_out}",sep="\t", quote = F)
#write.table(t,"${tvalue_out}",sep="\t", quote = F)
pdf("${heatmap_out}")
#pm <- pheatmap(cor,cluster_rows=T,cluster_cols=T,display_numbers=T,fontsize_number=4,number_format = "%.2f") #default euclidean
pm <- pheatmap(cor,cluster_rows=T,cluster_cols=T,display_numbers=T,fontsize_number=4) #default euclidean
tre_col <- as.phylo(pm$tree_col)
tre_row <- as.phylo(pm$tree_row)
write.tree(tre_col, "${col_tree}")
write.tree(tre_row, "${row_tree}")
dev.off()
