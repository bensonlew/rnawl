# load package

library(getopt)

# set options
command <- matrix(c(
    "json", "i", 1, "character", "setting file in JSON format",
    "help", "h", 0, "logical", "show this help message and exit"
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

if (!is.null(opts$help)) {
    cat(getopt(command, usage = TRUE))
    q(status = 2)
}
if (is.null(opts$json)) {
    cat(getopt(command, usage = TRUE))
    q(status = -1)
}

library(rjson)
library(DESeq2)

print("start loading required data")
json.lst <- fromJSON(file = opts$json)
count.path <- json.lst$count
sample.vec <- json.lst$samples
BATCH = json.lst$Batch
PAIRED = json.lst$Paired
paired.vec = json.lst$paireds
group.vec <- json.lst$groups
batch.ves = json.lst$batchs
contrast.vec <- json.lst$contrasts
together <- json.lst$together
output.path <- json.lst$output
count.df <- read.delim(count.path, row.names = 1)

if (together) {
    print("start computing differential expression")
    if (BATCH){
    col.df <- data.frame(group = group.vec,batch = batch.ves)
    rownames(col.df) <- sample.vec
    print(col.df)
    dds.i.obj <- DESeqDataSetFromMatrix(countData = count.df, colData = col.df, design = ~group+batch)
    } else if (PAIRED){
    col.df <- data.frame(group = group.vec,paired = paired.ves)
    rownames(col.df) <- sample.vec
    print(col.df)
    dds.i.obj <- DESeqDataSetFromMatrix(countData = count.df, colData = col.df, design = ~group+paired)
    } else if (PAIRED == FALSE & BATCH == FALSE) {
    col.df <- data.frame(group = group.vec)
    rownames(col.df) <- sample.vec
    print(col.df)
    dds.i.obj <- DESeqDataSetFromMatrix(countData = count.df, colData = col.df, design = ~group)
    }
    dds.o.obj <- DESeq(dds.i.obj)
    sf.vec <- dds.o.obj$sizeFactor
    count.norm.mat <- counts(dds.o.obj, normalized = TRUE)
    dr.obj.lst <- list()
    for (i in seq(length(contrast.vec))) {
        case <- contrast.vec[[i]]$case
        ctrl <- contrast.vec[[i]]$ctrl
        cmp <- paste(ctrl, case, sep = "_vs_")
        dr.obj.lst[[cmp]] <- results(dds.o.obj, contrast = c("group", case, ctrl))
    }
} else {
    count.df.lst <- list()
    col.df.lst <- list()
    sf.vec.lst <- list()
    count.norm.mat.lst <- list()
    dr.obj.lst <- list()
    for (i in seq(length(contrast.vec))) {
        case <- contrast.vec[[i]]$case
        ctrl <- contrast.vec[[i]]$ctrl
        cmp <- paste(ctrl, case, sep = "_vs_")
        print(paste("start computing differential expression at", cmp))
        case.idx.vec <- which(group.vec == case)
        ctrl.idx.vec <- which(group.vec == ctrl)
        case.group.vec <- group.vec[case.idx.vec]
        ctrl.group.vec <- group.vec[ctrl.idx.vec]
        case.sample.vec <- sample.vec[case.idx.vec]
        ctrl.sample.vec <- sample.vec[ctrl.idx.vec]
        col.df <- data.frame(group = c(case.group.vec, ctrl.group.vec))
        rownames(col.df) <- c(case.sample.vec, ctrl.sample.vec)
        col.df.lst[[cmp]] <- col.df
        count.df.lst[[cmp]] <- count.df[, c(case.sample.vec, ctrl.sample.vec)]
        dds.i.obj <- DESeqDataSetFromMatrix(countData = count.df.lst[[cmp]], colData = col.df.lst[[cmp]],
                                            design = ~group)
        dds.o.obj <- DESeq(dds.i.obj)
        sf.vec.lst[[cmp]] <- dds.o.obj$sizeFactor
        count.norm.mat <- counts(dds.o.obj, normalized = TRUE)
        count.norm.mat.lst[[cmp]] <- count.norm.mat
        dr.obj.lst[[cmp]] <- results(dds.o.obj, contrast = c("group", case, ctrl))
    }
}

print("start exporting results")
if (together) {
    write.table(sf.vec, file = file.path(output.path, "size_factor.txt"), quote = FALSE, sep = "\t", col.names = FALSE)
    write.table(t(c("seq_id", colnames(count.norm.mat))), file = file.path(output.path, "normalized.txt"),
                quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    write.table(count.norm.mat, file = file.path(output.path, "normalized.txt"), append = TRUE, quote = FALSE,
                sep = "\t", col.names = FALSE)
} else {
    for (cmp in names(sf.vec.lst)) {
        sf.vec <- sf.vec.lst[[cmp]]
        count.norm.mat <- count.norm.mat.lst[[cmp]]
        s.o.path <- file.path(output.path, paste(cmp, "size_factor", "txt", sep = "."))
        n.o.path <- file.path(output.path, paste(cmp, "normalized", "txt", sep = "."))
        write.table(sf.vec, file = s.o.path, quote = FALSE, sep = "\t", col.names = FALSE)
        write.table(t(c("seq_id", colnames(count.norm.mat))), file = n.o.path, quote = FALSE, sep = "\t",
                    row.names = FALSE, col.names = FALSE)
        write.table(count.norm.mat, file = n.o.path, append = TRUE, quote = FALSE, sep = "\t", col.names = FALSE)
    }
}

for (cmp in names(dr.obj.lst)) {
    o.path <- file.path(output.path, paste(cmp, "deseq2", "txt", sep = "."))
    write.table(t(c("seq_id", colnames(dr.obj.lst[[cmp]]))), file = o.path, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)
    write.table(dr.obj.lst[[cmp]], file = o.path, append = TRUE, quote = FALSE, sep = "\t", col.names = FALSE)
    print(paste("succeed in exporting", o.path))
}
