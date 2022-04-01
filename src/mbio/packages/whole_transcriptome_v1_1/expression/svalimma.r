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
library(sva)

# set variables

print("start loading required data")
json.lst <- fromJSON(file = opts$json)
count.path <- json.lst$count
group <- json.lst$groups
contrast.vec <- json.lst$contrasts
output.path <- json.lst$output
count.df <- read.delim(count.path, row.names = 1)
all_counts_filter <- count.df[which(rowSums(count.df)>0),]
all_counts_filter <- as.matrix(all_counts_filter)
mod <- model.matrix(~0+group)
mod0 <- cbind(rep(1,length(colnames(count.df))))
svobj <- svaseq(all_counts_filter, mod, mod0)
sv <- svobj$sv
sv <- as.matrix(sv)
colnames(sv) <- colnames(sv,do.NULL=FALSE,prefix = "sv_")
modSv <- cbind(mod,sv)

print("start getting DGEList and calculating norm factors")
d.obj <- DGEList(counts = count.df, group = group)
i.y.obj <- calcNormFactors(d.obj)
norm.count.df <- t(t(count.df) / i.y.obj$samples$norm.factors)

print("start exporting files containning norm information")
norm.factor.path <- file.path(output.path, "norm.factor.txt")
write.table(t(c("sample", "group", "lib.size", "norm.factors")), file = norm.factor.path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(i.y.obj$samples, file = norm.factor.path, append = TRUE, quote = FALSE, sep = "\t", col.names = FALSE)
norm.count.path <- file.path(output.path, "norm.count.txt")
write.table(t(c("seq_id", colnames(norm.count.df))), file = norm.count.path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(norm.count.df, file = norm.count.path, append = TRUE, quote = FALSE, sep = "\t", col.names = FALSE)


y_voom <- voom(i.y.obj, modSv)
fit <- lmFit(y_voom, modSv)
for (i in seq(length(contrast.vec))) {
    case <- contrast.vec[[i]]$case
    ctrl <- contrast.vec[[i]]$ctrl
    con = makeContrasts(paste(paste("group", case, sep=''), paste("group", ctrl, sep=''), sep = "-"), levels=modSv)
    fit2 <- contrasts.fit(fit, con)
    fit3 <- eBayes(fit2)
    table.obj <- topTable(fit3, n = nrow(count.df))
    o.path <- file.path(output.path, paste(paste(ctrl, case, sep = "_vs_"), "svaseqlimma", "txt", sep = "."))
    write.table(t(c("seq_id", colnames(table.obj))), file = o.path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    write.table(table.obj, file = o.path, append = TRUE, quote = FALSE, sep = "\t", col.names = FALSE)
    print(paste("succeed in exporting", o.path))
}

