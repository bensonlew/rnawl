options(warn=-100)
method<-'${method}'
k<- ${sub_num}
distance_method <- "${distance_method}"
input_matrix<-"${input_matrix}"
lognorm<-${lognorm}
cltype<-"${cltype}" #### both row column none

library(cluster)
library(Biobase)
library(mclust)
library(fpc)
library(ape)

data = read.delim(input_matrix, header=T, check.names=F, sep="\t")
rownames(data) = data[,1] # set rownames to gene identifiers
data = data[,2:length(data[1,])] # remove the gene column since its now the rowname value
data = as.matrix(data) # convert to matrix
colnames(data)<-colnames(data)

if(lognorm!=0){
    data = log(data+1,base=lognorm)
    final_data = t(scale(t(data))) # center and scale rows
    hc_genes = agnes(final_data, diss=FALSE, method = distance_method, metric='euclidean') # method: default-complete-linkage metric: default-euclidean cluster genes
    hc_samples = hclust(as.dist(1-cor(final_data, method="spearman")), method=distance_method) # cluster conditions
}
if(lognorm==0){
    final_data = t(scale(t(data)))
    hc_genes = agnes(final_data,diss=FALSE, method=distance_method, metric='euclidean') # cluster genes
    hc_samples = hclust(as.dist(1-cor(final_data, method="spearman")), method=distance_method) # cluster conditions
}

# if(k==0){
#     pamk <- pamk(final_data,diss=FALSE)
#     k <- pamk$nc
# }

####hclust
if ((method == "hclust") || (method == "both")){
    dir.create("hclust")

    ####### output the odered matrix after clustered
    if(cltype=="both"){
        order_mat<-data[hc_genes$order[nrow(data):1],hc_samples$order]
        write.table(order_mat,"hclust/hclust_heatmap.xls",sep="\t",col.names=T,row.names=T,quote=F)
        write.table(colnames(data)[hc_samples$order],"hc_sample_order",col.names=F,row.names=F,quote=F)
        write.table(rownames(data)[hc_genes$order[nrow(data):1]],"hc_gene_order",col.names=F,row.names=F,quote=F)
    }
    if(cltype=="row"){
        order_mat<-data[hc_genes$order[nrow(data):1],]
        write.table(order_mat,"hclust/hclust_heatmap.xls",sep="\t",col.names=T,row.names=T,quote=F)
        write.table(rownames(data)[hc_genes$order[nrow(data):1]],"hc_gene_order",col.names=F,row.names=F,quote=F)
     }
    if(cltype=="col"){
        order_mat<-data[,hc_samples$order]
        write.table(order_mat,"hclust/hclust_heatmap.xls",sep="\t",col.names=T,row.names=T,quote=F)
        write.table(colnames(data)[hc_samples$order],"hc_sample_order",col.names=F,row.names=F,quote=F)
     }

    ####### newick tree
    genes_tree <- as.phylo(as.hclust(hc_genes))
    samples_tree <- as.phylo(hc_samples)
    write.tree(genes_tree,"hclust/genes_tree.txt")
    write.tree(samples_tree,"hclust/samples_tree.txt")

    ####### the subclusters output
    if (k != 0){
        gene_partition_assignments <- cutree(as.hclust(hc_genes), k);
        # dir.create("hclust_subcluster")
        gene_names = rownames(final_data)
        num_cols = length(final_data[1,])
        for (i in 1:k) {
            partition_i = (gene_partition_assignments == i)
            partition_centered_data = final_data[partition_i,]
            # if the partition involves only one row, then it returns a vector instead of a table
            if (sum(partition_i) == 1){
                dim(partition_centered_data) = c(1,num_cols)
                colnames(partition_centered_data) = colnames(final_data)
                rownames(partition_centered_data) = gene_names[partition_i]
            }
            outfile = paste("hclust", "/subcluster_", i, sep='')
            write.table(partition_centered_data, file=outfile, quote=F, sep="\t")
        }
    }
}

####kmeans
if ((method == "kmeans") || (method == "both")){
    dir.create("kmeans")

    ###plots: within groups sum of squares
    if (k != 0){
        km <- kmeans(final_data, k)
        m1 <- cbind(km$cluster,final_data)
        m2 <- order(m1[, 1])
        m3 <- m1[m2, ]
        m4 <- m3[,-1]
        ###kmeans clusters
        write.table(m3,"kmeans/kmeans_heatmap.xls",sep="\t",col.names=T,row.names=T,quote=F)

        ###kmeans subclusters
        # dir.create("kmeans_subcluster")
            gene_names = rownames(final_data)
            num_cols = length(final_data[1,])
        for (i in 1:k) {
            sub_i = (km$cluster == i)
            sub_data = final_data[sub_i,]
            if (sum(sub_i) == 1){
                dim(sub_data) = c(1,num_cols)
                colnames(sub_data) = colnames(final_data)
                rownames(sub_data) = gene_names[sub_i]
            }
            outfile = paste("kmeans", "/subcluster_", i, sep='')
            write.table(sub_data, file=outfile, quote=F, sep="\t")
        }
    }
}
