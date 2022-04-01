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
library(edgeR)

# set variables

print("start loading required data")
json.lst <- fromJSON(file = opts$json)
count.path <- json.lst$count
group.vec <- json.lst$groups
BATCH = json.lst$Batch
PAIRED = json.lst$Paired
batch.vec = json.lst$batchs
paired.vec = json.lst$paireds
contrast.vec <- json.lst$contrasts
output.path <- json.lst$output
fit.method <- json.lst$fit_method
count.df <- read.delim(count.path, row.names = 1)
count.df = floor(count.df)

print("start creating group.factor")
group = as.factor(group.vec)

if (BATCH == TRUE) {
print("start creating batch.factor")
batch = as.factor(batch.vec)
}

if (PAIRED == TRUE) {

print("start creating paired.factor")
paired = as.factor(paired.vec)
}


print("start getting DGEList and calculating norm factors")
d.obj <- DGEList(counts = count.df, group = group)
i.y.obj <- calcNormFactors(d.obj)
norm.count.df <-t(t(count.df) / i.y.obj$samples$norm.factors)

print("start exporting files containning norm information")
norm.factor.path <- file.path(output.path, "norm.factor.txt")
write.table(t(c("sample", "group", "lib.size", "norm.factors")), file = norm.factor.path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(i.y.obj$samples, file = norm.factor.path, append = TRUE, quote = FALSE, sep = "\t", col.names = FALSE)
norm.count.path <- file.path(output.path, "norm.count.txt")
write.table(t(c("seq_id", colnames(norm.count.df))), file = norm.count.path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(norm.count.df, file = norm.count.path, append = TRUE, quote = FALSE, sep = "\t", col.names = FALSE)

print("start creating design matrix")
if (BATCH == FALSE & PAIRED == FALSE) {
design.mat = model.matrix(~0+group)
design.mat
#colnames(design.mat)=levels(group)
#rownames(design.mat) <- colnames(count.df)
} else if (BATCH == TRUE & PAIRED == FALSE) {
design.mat = model.matrix(~0+group+batch)
design.mat
#colnames(design.mat)=levels(group)
#rownames(design.mat) <- colnames(count.df)
} else if (BATCH == FALSE & PAIRED == TRUE) {
design.mat = model.matrix(~0+group+paired)
design.mat
}
#colnames(design.mat) <- levels(i.y.obj$samples$group)

print("start estimating dispersions across all tags")
o.y.obj <- estimateDisp(i.y.obj, design.mat)

print("start fitting NB glm model to count data")
if (fit.method == "glmFit") {
    fit.obj <- glmFit(o.y.obj, design.mat)
} else if (fit.method == "glmQLFit") {
    fit.obj <- glmQLFit(o.y.obj, design.mat)
}

print("start conducting genewise statistical tests for given contrasts")
for (i in seq(length(contrast.vec))) {
    case = contrast.vec[[i]]$case
    ctrl = contrast.vec[[i]]$ctrl
    group_case <- paste('group',contrast.vec[[i]]$case,sep='')
    group_ctrl <- paste('group',contrast.vec[[i]]$ctrl,sep='')

    if (fit.method == "exactTest") {
        qlf.obj <- exactTest(o.y.obj, pair = c(ctrl, case), dispersion = 0.1)
    } else if (fit.method == "glmFit") {
        qlf.obj <- glmLRT(fit.obj, contrast = makeContrasts(paste(group_case, group_ctrl, sep = "-"), levels = design.mat))
    } else if (fit.method == "glmQLFit") {
        qlf.obj <- glmQLFTest(fit.obj, contrast = makeContrasts(paste(group_case, group_ctrl, sep = "-"), levels = design.mat))
    }
    table.obj <- topTags(qlf.obj, n = dim(count.df)[1])
    o.path <- file.path(output.path, paste(paste(ctrl, case, sep = "_vs_"), "edger", "txt", sep = "."))
    write.table(t(c("seq_id", colnames(table.obj))), file = o.path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    write.table(table.obj, file = o.path, append = TRUE, quote = FALSE, sep = "\t", col.names = FALSE)
    print(paste("succeed in exporting", o.path))
}
