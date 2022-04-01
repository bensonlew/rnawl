#!/usr/bin/env Rscript
times<-Sys.time();
library('getopt')
library(pheatmap)
options(bitmapType='cairo')
spec = matrix(c(
	'i','a',0,'character',
	'o','b',0,'character',
	'help','c', 0, 'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	
Usage:
	--i different matrix between samples
	--o figure name
	--help		usage
\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$i) ) { print_usage(spec) }
if ( is.null(opt$o) ) { print_usage(spec) }



dat_raw <- read.table(opt$i, sep = "\t", header = T,check.names=F)
rownames(dat_raw)<-dat_raw[,1]
dat_raw<-dat_raw[,-1]
dat_raw<-t(dat_raw)
outpdf<-paste(opt$o,".pdf",sep="")
pheatmap(dat_raw,cluster_rows=T,cluster_cols=T, 
	fontsize_number = 8,	
	filename=outpdf,
)
outpng<-paste(opt$o,".png",sep="")
pheatmap(dat_raw,cluster_rows=T,cluster_cols=T, 
	fontsize_number = 8,
	filename=outpng,
)
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
