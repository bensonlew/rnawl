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
library(DEGseq)

# set variables

print("start loading required data")
json.lst <- fromJSON(file = opts$json)
count.path <- json.lst$count
contrast.lst <- json.lst$contrasts
column.lst <- json.lst$columns
output.path <- json.lst$output

print("start computing differential expression")
for (i in seq(length(column.lst))) {
    case.col.vec <- column.lst[[i]]$case
    ctrl.col.vec <- column.lst[[i]]$ctrl
    case.label <- contrast.lst[[i]]$case
    ctrl.label <- contrast.lst[[i]]$ctrl
    case.exp.mat <- readGeneExp(file = count.path, geneCol = 1, valCol = case.col.vec)
    ctrl.exp.mat <- readGeneExp(file = count.path, geneCol = 1, valCol = ctrl.col.vec)
    output.dir <- file.path(output.path, paste(ctrl.label, case.label, sep = "_vs_"))
    DEGexp(geneExpMatrix1 = case.exp.mat, geneCol1 = 1, expCol1 = seq(dim(case.exp.mat)[2], from = 2), groupLabel1 = case.label,
           geneExpMatrix2 = ctrl.exp.mat, geneCol2 = 1, expCol2 = seq(dim(ctrl.exp.mat)[2], from = 2), groupLabel2 = ctrl.label,
           method = "MARS", outputDir = output.dir)
    print(paste("succeed in exporting results to", output.dir))
}
