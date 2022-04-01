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
library(GSVA)
library(GSEABase)

json.lst <- fromJSON(file = opts$json)
matrix.path <-json.lst$matrix
gmx <- json.lst$gmx
method <- json.lst$method
kcdf <- json.lst$kcdf
min_num <- json.lst$min_num
max_num <- json.lst$max_num
es_score <- json.lst$es_score
output <- json.lst$output

geneset_list <- getGmt(gmx)
matrix <- read.table(matrix.path, row.names = 1, header = TRUE, sep='\t')
matrix <- as.matrix(matrix)
gsva_es <- gsva(matrix, geneset_list, method = method, min.sz = min_num, max.sz = max_num)
if (method=='ssgsea'){
gsva_es = gsva_es
}else{
gsva_es = gsva_es$es.obs
}
col = c('geneset' ,colnames(gsva_es))
write.table(t(col), file = output, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(gsva_es, file = output,append = TRUE, quote = FALSE, sep = "\t", col.names = FALSE, row.names = TRUE)




