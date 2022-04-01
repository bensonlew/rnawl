#!/usr/bin/R
#fpkm="/mnt/ilustre/users/sanger-dev/workspace/20170602/GenesetVenn_tsg_2000_9965_3754/geneset_file_geneset_venn"
#outfile_path = "/mnt/ilustre/users/sanger-dev/workspace/20170602/GenesetVenn_tsg_2000_9965_3754/GenesetVenn/output"
fpkm = "/mnt/ilustre/users/sanger-dev/workspace/20170602/GenesetVenn_tsg_2000_2210_9365/geneset_file_geneset_venn"
outfile_path = "/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/ref_rna/express"
data = read.delim(fpkm, sep="\t",header=F,as.is=T)
lst = list()
for(i in 1:dim(data)[1]){
    seq_id = strsplit(data[i,2],split=",")
    lst[[data[i,1]]] = as.vector(seq_id[[1]])
}

intersect_n  <- function(x){
    n <-length(x)
    x2 <- x[[1]]
    for (i in 2: n){
        x2 <- intersect(x[[i]], x2)
    }
    return (x2)
}

getset <- function(lst){
    # lst is a list
    sn <-length(names(lst))
    sets <- list()
    sls <- list()
    maxl <-0
    # get all intersect
    for (i in sn: 2){
        sl <- combn(names(lst), i, simplify=FALSE)
        inter <-lapply(sl, function(x)intersect_n(lst[x]))
        names(inter) <-lapply(sl, function(x)paste(x, collapse=" & "))
        sets <- c(sets, inter)
        sls <- c(sls, sl)
    }
    # get all unique
    for (i in 1: sn){
        uniq <- list(setdiff(unlist(lst[names(lst)[i]]), unlist(lst[names(lst)[-i]])))
        names(uniq) <- paste(names(lst)[i]," only")
        sets <- c(sets, uniq)
    }
    return (sets)
}

sets <- getset(lst)
print(names(sets))

lst_len <- length(sets)
venn_matrix <- matrix(, lst_len, 3)
for(i in 1:lst_len){
    venn_matrix[i,1] <- names(sets[i])
    venn_matrix[i,2] <- length(sets[[i]])
    venn_matrix[i,3] <- paste(sets[[i]],collapse =",")
}
otuset <- paste(outfile_path,"/venn_table.xls",collapse="",sep="")
write.table(venn_matrix,otuset,sep = "\t",eol="\n",row.names=FALSE,col.names=FALSE,quote=FALSE)
