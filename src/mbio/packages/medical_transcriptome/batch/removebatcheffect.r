# load package

library(getopt)

# set options
command <- matrix(c(
    "json", "i", 1, "character", "setting file in JSON format",
    "help", "h", 0, "logical", "show this help message and exit"
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

if (! is.null(opts$help)) {
    cat(getopt(command, usage = TRUE))
    q(status = 2)
}
if (is.null(opts$json)) {
    cat(getopt(command, usage = TRUE))
    q(status = - 1)
}

library(rjson)
library(limma)

# set variables

print("start loading required data")
json.lst <- fromJSON(file = opts$json)
count.path <- json.lst$count
group.vec <- json.lst$groups
batch.vec = json.lst$batchs
output.path <- json.lst$output
count.df = read.delim(count.path, row.name=1)
#count_df = count.df[which(rowSums(count.df)>0),]
batch = as.numeric(batch.vec)
group = as.factor(group.vec)
pheno = cbind(batch, group)
pheno = as.data.frame(pheno)
design = matrix(pheno$group,length(pheno$group),1)
print("start batch effects analysis")
col_name = colnames(count.df)
write.table(t(c('seq_id',col_name)), file = file.path(output.path, 'count_batch.txt'), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
exp_batch = removeBatchEffect(count.df, batch = pheno$batch, design = design)
exp_batch[exp_batch<0]=0
write.table(exp_batch, file = file.path(output.path, 'count_batch.txt'), sep = '\t', quote = FALSE, row.names = TRUE, col.name = FALSE, append = TRUE)


