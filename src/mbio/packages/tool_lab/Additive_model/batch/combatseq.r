# load package

library(getopt)
source('/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/tool_lab/Additive_model/combat_seq.R')
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
#library(sva)


# set variables

print("start loading required data")
json.lst <- fromJSON(file = opts$json)
count.path <- json.lst$count
group.vec <- json.lst$groups
output.path <- json.lst$output
batch = json.lst$batchs
count.df <- read.delim(count.path, row.names = 1)
count.df = floor(count.df)
count.df = as.matrix(count.df)

adjusted_counts <- ComBat_seq(count.df, batch=batch, group=group.vec)

col_name = colnames(count.df)
write.table(t(c('seq_id',col_name)), file = file.path(output.path, 'count_combat_seq.txt'), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(adjusted_counts, file = file.path(output.path, 'count_combat_seq.txt'), sep = '\t', quote = FALSE, row.names = TRUE, col.name = FALSE, append = TRUE)


