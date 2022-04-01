#######################################################################################################################################转换了蛋白质ID
library(STRINGdb)
library(dplyr)
library(igraph)
args <- commandArgs(T)
#if(length(args) < 5)
#{
#            stop("Rscript  PPInetwork_predict.r <gene_list.txt> <species> <outdir> <combined_score> <score_cutoff>")
#}
##四个参数，args[1]输入map后数据；args[2]物种名称如9606;args[3]输出路径;args[4]combin_score值(按照combined_score进行降序排序，然后取前300行进行计算)
cat(args[1])
cat("\n")
cat(args[2])
spec <- as.numeric(args[2])
cat("\n")
outdir <- args[3]
cat(args[3])
dir.create(outdir)
cat("\n")
combinescore <- as.numeric(args[4])
cat(args[4])
cat("\n")
database <- args[5]
score_cutoff <- as.numeric(args[6])
cat(args[6])
cat("\n")
score_cutoff <- as.double(args[6]) * 1000

string_db<-STRINGdb$new(version="10", species = spec, input_directory=database)
diff_exp_example1 <- read.table(args[1], header = T, sep = "\t")
diff_exp_example1 <- diff_exp_example1[!duplicated(diff_exp_example1["STRING_id"]), ]
inter <- string_db$get_interactions(diff_exp_example1$STRING_id)
#######################输出所有的互作组
protein <- string_db$get_proteins()
write.table(protein, file = paste(outdir,"/gene_protein.txt", sep = ""),sep = "\t", row.names = FALSE, quote = FALSE) #输出蛋白质信息
protein_12 <- protein[1:2]
#获得蛋白与logFC对应文件
#######################输出combine_score>args[4] 的互作组
inter_x <- unique(inter)
inter_x12 <- cbind( inter_x[1:15], inter_x[16])
merg1_x <- merge(inter_x12, protein_12, by.x = "from", by.y = "protein_external_id", all.x = TRUE)
merg_x <- merge(merg1_x, protein_12, by.x = "to", by.y = "protein_external_id", all.x = TRUE)
x_na = is.na(merg_x$preferred_name.x) | merg_x$preferred_name.x == ""
y_na = is.na(merg_x$preferred_name.y) | merg_x$preferred_name.y == ""
merg_x$preferred_name.y[y_na] = merg_x$to[y_na]
merg_x$preferred_name.x[x_na] = merg_x$from[x_na]
PPI_data_x <- select(merg_x, from = preferred_name.x, to = preferred_name.y, combined_score = combined_score)
sorted_ppi_data <- PPI_data_x[order(PPI_data_x[,3], decreasing = T),]   #按照combined_score进行降序排序
sorted_merg_x <- merg_x[order(merg_x[,16], decreasing = T),]   #按照combined_score进行降序排序
sorted_ppi_data = sorted_ppi_data[sorted_ppi_data$combined_score > score_cutoff ,]
sorted_merg_x = sorted_merg_x[sorted_merg_x$combined_score > score_cutoff ,]
data_lines <- nrow(sorted_ppi_data)
if(combinescore > data_lines){n <- data_lines}else if(combinescore <= data_lines){n <- combinescore}
PPI_data_x <- sorted_ppi_data[0:n,]
sorted_merg_x <- sorted_merg_x[0:n, ]
PPI_data_x <- select(sorted_merg_x, from = preferred_name.x, to = preferred_name.y, combined_score = combined_score)
sorted_merg_x[3:16]=sorted_merg_x[3:16]/1000
#if (nrow(PPI_data_x) < 1) {stop("不存在网络的边文件，无法进行网络分析！")}
#######################计算网络的测度
myarray1 <- array(PPI_data_x$from)
myarray2 <- array(PPI_data_x$to)
myarray <- append(myarray1, myarray2)
myarray_sort <- sort(myarray)
myarray_sort_uniq <- unique(myarray_sort)
myarray_sort_uniq <- myarray_sort_uniq[myarray_sort_uniq != ""]
actors <- data.frame(name = myarray_sort_uniq)
relations <- data.frame(from = PPI_data_x$from, to = PPI_data_x$to)
g <- graph.data.frame(relations, directed = F, vertices = actors)
all_nodes1 <- data.frame(node = myarray_sort_uniq, degree = degree(g))
all_nodes2 <- merge(all_nodes1, protein_12, by.x = "node", by.y = "preferred_name", all.x = TRUE)
ex_na = is.na(all_nodes2$protein_external_id)
all_nodes2$protein_external_id[ex_na]=as.vector(all_nodes2$node[ex_na])
#all_nodes3 <- merge(all_nodes2, diff_exp_example1, by.x = "protein_external_id", by.y = "STRING_id", all.x = TRUE)
all_nodes3 <- merge(all_nodes2, diff_exp_example1, by.x = "protein_external_id", by.y = "STRING_id")
all_nodes <- select(all_nodes3, node = node, degree = degree, accession_id = accession_id, STRING_id = protein_external_id, method=method)
#######################统计网络的特征值
num_of_nodes <- length(V(g))
num_of_edges <- length(E(g))
cluster_coefficient <- transitivity(g, type = "average")
transitivity <- transitivity(g)
average_node_degree <- sum(all_nodes$degree)/num_of_nodes
average_path_length <- average.path.length(g)
network_stats <- rbind(num_of_nodes,num_of_edges,average_node_degree,average_path_length,cluster_coefficient,transitivity)
write.table(network_stats, file = paste(outdir,"/network_stats.txt", sep = ""), sep = "\t", col.name = FALSE, quote = FALSE)
acc2string=data.frame(diff_exp_example1$accession_id)
rownames(acc2string)=diff_exp_example1$STRING_id
sorted_merg_x$from_accession = acc2string[sorted_merg_x$from,]
sorted_merg_x$to_accession = acc2string[sorted_merg_x$to,]
write.table(sorted_merg_x, file = paste(outdir,"/interaction_detail.txt", sep = ""), sep = "\t", row.names=FALSE, col.name = TRUE, quote = FALSE)
protein=protein[c(1,3,4)]
colnames(protein) = c('STRING_id', 'STRING_length', 'STRING_description')
all_nodes=merge(all_nodes, protein, by='STRING_id', all.x=TRUE)
write.table(all_nodes, file = paste(outdir,"/all_nodes.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(PPI_data_x, file = paste(outdir,"/interaction.txt", sep = ""),sep = "\t", row.names = FALSE, quote = FALSE)
