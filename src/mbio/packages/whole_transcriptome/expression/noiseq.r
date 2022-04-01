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
library(NOISeq)

# set variables

print("start loading required data")
json.lst <- fromJSON(file = opts$json)
count.path <- json.lst$count
group <- json.lst$groups
contrast.vec <- json.lst$contrasts
singal_sample <- json.lst$singal_sample
output.path <- json.lst$output
count.df <- read.delim(count.path, row.names = 1)

tmp_col <- colnames(count.df)
myfactors <- data.frame(group, row.names=tmp_col)

print("start exporting files containning norm information")
norm.count.path <- file.path(output.path, "norm.count.txt")
normalized_count <- tmm(count.df, long = 1000, lc = 0, k = 0.5)
colnames <- c("seq_id",colnames(normalized_count))
write.table(t(colnames), file=norm.count.path , sep="\t", quote=F, col.names=F,row.names=F)
write.table(normalized_count, file=norm.count.path, sep="\t", quote=F, row.names=T, col.names=F, append=T)

print("start calculate")
mydata <- readData(data = normalized_count, factors = myfactors)
for (i in seq(length(contrast.vec))) {
     case <- contrast.vec[[i]]$case
     ctrl <- contrast.vec[[i]]$ctrl
     if (singal_sample == 'no') {
        mynoiseq <- noiseqbio(mydata, norm = "n", factor = "group", conditions = c(case, ctrl))
     }else if (singal_sample == 'yes') {
        mynoiseq <- noiseq(mydata, norm = "n", replicates = "no", factor = "{}", conditions = c("{}", "{}"))
     }
     results <- mynoiseq@results[[1]]

     colnames <- c("seq_id",colnames(results))
     o.path <- file.path(output.path, paste(paste(ctrl, case, sep = "_vs_"), "NOIseq", "txt", sep = "."))
     write.table(t(colnames), file=o.path, sep="\t", quote=F, col.names=F,row.names=F)
     write.table(results, file=o.path, sep="\t", quote=F, row.names=T, col.names=F, append=T)
}





















