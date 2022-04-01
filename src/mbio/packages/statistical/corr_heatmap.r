library(pheatmap)
library(psych)
library(ape)
data=read.table("${inputfile}",header=T,row.names=1,comment.char = '',sep="\t")
#pmcorr <- pheatmap(data,cluster_rows=T,cluster_cols=T,display_numbers=T,fontsize_number=4,number_format = "%.2f")

#col_pmcorr <- pheatmap(data,clustering_method="${col_cluster_method}",silent=TRUE,cluster_rows=T,cluster_cols=T,display_numbers=T,fontsize_number=4,number_format = "%.2f")
#row_pmcorr <- pheatmap(data,clustering_method="${row_cluster_method}",silent=TRUE,cluster_rows=T,cluster_cols=T,display_numbers=T,fontsize_number=4,number_format = "%.2f")
col_cluster_method <- "${col_cluster_method}"
if(col_cluster_method != ""){
    col_pmcorr <- pheatmap(data,clustering_method=col_cluster_method,silent=TRUE,cluster_rows=T,cluster_cols=T,display_numbers=T,fontsize_number=4)
    corr_tre_col <- as.phylo(col_pmcorr$tree_col)
    write.tree(corr_tre_col, "${col_tree}")
}
row_cluster_method <- "${row_cluster_method}"
if(row_cluster_method != ""){
    row_pmcorr <- pheatmap(data,clustering_method=row_cluster_method,silent=TRUE,cluster_rows=T,cluster_cols=T,display_numbers=T,fontsize_number=4)
    corr_tre_row <- as.phylo(row_pmcorr$tree_row)
    write.tree(corr_tre_row, "${row_tree}")
}
