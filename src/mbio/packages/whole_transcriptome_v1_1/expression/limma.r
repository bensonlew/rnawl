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
library(limma)

# set variables

print("start loading required data")
json.lst <- fromJSON(file = opts$json)
count.path <- json.lst$count
group.vec <- json.lst$groups
contrast.vec <- json.lst$contrasts
is_batch <- json.lst$is_batch
if (is_batch == TRUE){
    batch.vec <- json.lst$batchs
    batch_factor <- as.factor(batch_factor)
}
singal_sample <- json.lst$singal_sample
output.path <- json.lst$output
count.df <- read.delim(count.path, row.names = 1)

print("start getting DGEList and calculating norm factors")
d.obj <- DGEList(counts = count.df, group = group.vec)
i.y.obj <- calcNormFactors(d.obj)
norm.count.df <- t(t(count.df) / i.y.obj$samples$norm.factors)

print("start exporting files containning norm information")
norm.factor.path <- file.path(output.path, "norm.factor.txt")
write.table(t(c("sample", "group", "lib.size", "norm.factors")), file = norm.factor.path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(i.y.obj$samples, file = norm.factor.path, append = TRUE, quote = FALSE, sep = "\t", col.names = FALSE)
norm.count.path <- file.path(output.path, "norm.count.txt")
write.table(t(c("seq_id", colnames(norm.count.df))), file = norm.count.path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(norm.count.df, file = norm.count.path, append = TRUE, quote = FALSE, sep = "\t", col.names = FALSE)

print("if is a singal_sample")
if (singal_sample == 'no'){
    print("start creating design matrix")
    if (is_batch == True){
        design.mat <- model.matrix(~ 0 + group + batch_factor, data = i.y.obj$samples)
    }else{
        design.mat <- model.matrix(~ 0 + group, data = i.y.obj$samples)
    }

    y_voom <- voom(i.y.obj, design.mat)
    fit <- lmFit(y_voom, design.mat)
    for (i in seq(length(contrast.vec))) {
        case <- contrast.vec[[i]]$case
        ctrl <- contrast.vec[[i]]$ctrl
        con = makeContrasts(paste(paste("group", case, sep=''), paste("group", ctrl, sep=''), sep = "-"), levels=design.mat)
        fit2 <- contrasts.fit(fit, con)
        fit3 <- eBayes(fit2)
        table.obj <- topTable(fit3, n = nrow(count.df))
        o.path <- file.path(output.path, paste(paste(ctrl, case, sep = "_vs_"), "limma", "txt", sep = "."))
        write.table(t(c("seq_id", colnames(table.obj))), file = o.path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
        write.table(table.obj, file = o.path, append = TRUE, quote = FALSE, sep = "\t", col.names = FALSE)
        print(paste("succeed in exporting", o.path))
    }
}else if (singal_sample == "yes"){
        for (i in seq(length(contrast.vec))) {
        case <- contrast.vec[[i]]$case
        ctrl <- contrast.vec[[i]]$ctrl
        result <- exactTest(o.y.obj, pair = c(ctrl, case), dispersion = 0.1)
        table.obj <- topTable(result, n = nrow(count.df))
        o.path <- file.path(output.path, paste(paste(ctrl, case, sep = "_vs_"), "limma", "txt", sep = "."))
        write.table(t(c("seq_id", colnames(table.obj))), file = o.path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
        write.table(table.obj, file = o.path, append = TRUE, quote = FALSE, sep = "\t", col.names = FALSE)
        print(paste("succeed in exporting", o.path))
    }
}

