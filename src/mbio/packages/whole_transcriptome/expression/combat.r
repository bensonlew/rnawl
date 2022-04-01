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
exp.path <- json.lst$exp
pheno.path <- json.lst$pheno  #batch,condition
output.path <- json.lst$output
exp = read.table(exp.path, sep='\t', header=TRUE)
pheno = read.table(pheno.path, row.names=1, sep='\t', header=TRUE)

print("start batch effects analysis")
col_name = colnames(exp)
write.table(t(col_name), file = file.path(output.path, 'exp_batch.txt'), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
a = exp[,1]
row.names(exp) = a
exp = exp [,-1]
exp = as.matrix(exp)
model <- model.matrix(~condition, data = pheno)
combat_exp <- ComBat(dat = exp, batch = pheno$batch, mod = model)
write.table(combat_exp, file = file.path(output.path, 'exp_batch.txt'), sep = '\t', quote = FALSE, row.names = TRUE, col.name = FALSE, append = TRUE)


