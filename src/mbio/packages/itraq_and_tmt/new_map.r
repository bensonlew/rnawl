######################################################################################################################################
#将Ensemble ID map 到string数据库中的ID
library(STRINGdb)
args <- commandArgs(T)
#args[1]输入差异表达详情表；args[2]物种；args[3] FDR args[4]输出路径
#Rscript  map.r <diff_exp.xls> <species> <outdir>
cat(args[1])
cat("\n")
cat(args[2])
spec <- as.numeric(args[2])
cat("\n")
outdir <- args[3]
cat(args[3])
database<- args[4]
dir.create(outdir)
cat("\n")
string_db<-STRINGdb$new(version="10", species = spec, input_directory=database)
diff_exp_example1 <- read.table(args[1], header = T, sep = "\t")
#diff_exp_example1 <- diff_exp_example1[which(diff_exp_example1$FDR < FDR),]
data_mapped <- string_db$map(diff_exp_example1, "accession_id", takeFirst=FALSE, removeUnmappedRows=TRUE, quiet = TRUE)
#data_mapped <- string_db$map(diff_exp_example1, "seq_id", quiet = TRUE) #保留没有map上的差异基因
write.table(data_mapped, file = paste(outdir,"/seq_mapped.txt", sep = ""),sep = "\t", row.names = FALSE, quote = FALSE)
unmapped_gene <- setdiff(diff_exp_example1$accession_id, data_mapped$accession_id)
write.table(unmapped_gene, file = paste(outdir,"/unmapped_seq.txt", sep = ""),sep = "\t", row.names = FALSE, col.names=FALSE,quote = FALSE)
unmapped_db <- setdiff(string_db$proteins$protein_external_id, data_mapped$STRING_id)
write.table(unmapped_db, file = paste(outdir,"/unmapped_db.txt", sep = ""),sep = "\t", row.names = FALSE, col.names=FALSE,quote = FALSE)
