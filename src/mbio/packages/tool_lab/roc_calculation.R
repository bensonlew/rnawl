#!/usr/bin/env Rscript
# author: binbin.zhao@20200708
times<-Sys.time()
library('getopt');
library(pROC);
print("---------------------------");
spec = matrix(c(
	'in','i',1,'character',
	'od','o',1,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4);
print(spec);

opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("Usage example: \n")
	cat("
Usage example:
	Rscript indel_len.r --i  --o

Usage:
	--i     insert_size file
	--od	the output dir
	--help		usage
\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$i) ) { print_usage(spec) }
if ( is.null(opt$o) ) { print_usage(spec) }
data_roc<-read.table(opt$i,header=T)
attach(data_roc)
for(i in 1:(length(data_roc[0,])-1)){
	print(i)
	print("**********************")
	print(data_roc[,1])
	print("///////////////////////////")
	m <- roc(data_roc[,1], data_roc[,i+1], plot = T)
	data <- coords(m)
	file=paste(opt$o,"/", colnames(data_roc)[i+1],sep="")
	write.table(data, file)
	write.table(coords(m, "best"), file, append = TRUE)
	write.table(auc(m), file, append = TRUE)
}

detach(data_roc)
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)

