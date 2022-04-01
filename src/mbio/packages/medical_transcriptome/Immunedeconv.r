

## Reference https://icbi-lab.github.io/immunedeconv/articles/immunedeconv.html
#get parameters
library(getopt)

command = matrix(c(
    'help', 'h',  0, 'loical',
    'count', 'c', 1, 'character',
    'method', 'm', 1, 'character',
    'tumor', 't', 1, 'character',
    'cibersort', 'b', 1, 'character',
    'LM22', 'l', 1, 'character',
    'probesets', 'p', 1, 'character',
    'genes', 'g', 1, 'character',
    'out', 'o', 2,  'character'
  ),byrow=TRUE,ncol =4);
opt = getopt(command);
print_usage <- function(command=NULL){
        cat(getopt(command, usage=T));

        cat("
Usage example: 
Options:
        --help    NULL      get this help
        --count  <file>  a count matrixs of genes (rows) vs. samples (columns) .
        --method  <str>  default quantiseq, Alternative methods: 'quantiseq', 'timer',  'mcp_counter', 'xcell', 'epic'.
        --tumor  <str>  specify the tumor type of samples. if set method 'timer' , it can not ommit.
        --out  <dir>  the output dir, if not exist,it will be created.
        \n")
 q(status=1);
}
if ( is.null(opt$count) ) { print_usage(command) }
if ( is.null(opt$method) ){opt$method ="quantiseq" }
if ( !is.null(opt$tumor)){
    Tumor <- read.table(opt$tumor,sep="\t",header=F)
    Tumor <- unlist(Tumor[1,],use.names=F)}

##load library
library(dplyr)
library(ggplot2)
library(tidyr)
library(immunedeconv)
library(tibble)

Count <- read.table(opt$count,header=T,sep="\t",row.names=1)

if (opt$method == "cibersort"){
    set_cibersort_binary(opt$cibersort)
    set_cibersort_mat(opt$LM22)
}else if (opt$method == "timer"){
    Res <- deconvolute(Count,"timer",indications=Tumor)
}else if (opt$method == "mcp_counter"){
    Res <- deconvolute(Count, "mcp_counter", probesets=read.table(opt$probesets,sep="\t",stringsAsFactors=FALSE,colClasses="character"), genes=read.table(opt$genes,sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE))
}else {
    Res <- deconvolute(Count, opt$method)
}

setwd(opt$out)
Res <- as.data.frame(Res)
write.table(Res,paste0(opt$method,"_Immu_cell.txt"),row.names=F, quote=F, sep='\t')
