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
library(sva)

# set variables

print("start loading required data")
json.lst <- fromJSON(file = opts$json)
count.path <- json.lst$count
group.vec <- json.lst$groups
batch.vec = json.lst$batchs
output.path <- json.lst$output
counts = read.delim(count.path, row.name=1)
counts = as.matrix(counts)
group = as.factor(group.vec)
batch = as.factor(batch.vec)

## Remove genes with only 0 counts in any batch
keep_batch <- lapply(levels(batch), function(b){
  if (class(counts[,batch==b])=='numeric'){
  which(apply(matrix(counts[, batch==b]), 1, function(x){!all(x==0)}))
  } else{
  which(apply(counts[, batch==b], 1, function(x){!all(x==0)}))
}})
keep_batch <- Reduce(intersect, keep_batch)
keep_group = lapply(levels(group), function(b){
  if (class(counts[,group==b])=='numeric'){
  which(apply(matrix(counts[, group==b]), 1, function(x){!all(x==0)}))
  } else{
  which(apply(counts[, group==b], 1, function(x){!all(x==0)}))
}})
keep_group = Reduce(intersect, keep_group)
keep = intersect(keep_batch,keep_group)
rm <- setdiff(1:nrow(counts), keep)
countsOri <- counts
counts_batch_group <- counts[keep, ]


print("start batch effects analysis")
col_name = colnames(counts)
write.table(t(c('seq_id',(col_name))), file = file.path(output.path, 'count_batch.txt'), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

model <- model.matrix(~group)
combat_exp <- ComBat(dat = counts_batch_group, batch = batch, mod = model)

## Add back genes with only 0 counts in any batch (so that dimensions won't change)
adjust_counts_whole <- matrix(NA, nrow=nrow(countsOri), ncol=ncol(countsOri))
dimnames(adjust_counts_whole) <- dimnames(countsOri)
adjust_counts_whole[keep, ] <- combat_exp
adjust_counts_whole[rm, ] <- countsOri[rm, ]
adjust_counts_whole[adjust_counts_whole<0]=0
write.table(adjust_counts_whole, file = file.path(output.path, 'count_batch.txt'), sep = '\t', quote = FALSE, row.names = TRUE, col.name = FALSE, append = TRUE)


