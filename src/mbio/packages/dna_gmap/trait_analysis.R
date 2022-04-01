#!/usr/bin/env Rscript
# load library
times<-Sys.time()
library(ggplot2)
# library(ggpubr)
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'help' , 'h', 0, "logical",
	'infile' , 'i', 1, "character",
	'outfile' , 'o', 1, "character"
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
# define usage function
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example:
Options:
	--help		NULL 		get this help
	--infile 	character 	the input file [forced]
	--outfile 	character 	the filename for output graph [forced]
	\n")
	q(status=1);
}
if ( is.null(opt$infile) )	{ print_usage(spec) }
if ( is.null(opt$outfile) )	{ print_usage(spec) }

Trait <- read.table(opt$infile,header=T,comment.char="^")
TraitName=colnames(Trait)
data<-NULL;
for (i in 2:length(TraitName)){
	new<-ordered(na.omit(Trait[,i]));
	min<-min(new);
	max<-max(new);
	mean<-mean(new);
	var<-sd(new);
	num<-length(new);
	normt<-shapiro.test(as.numeric(new))
	# min<-min(na.omit(Trait[,i]));
	# max<-min(na.omit(Trait[,i]));
	# mean<-mean(na.omit(Trait[,i]));
	# var<-sd(na.omit(Trait[,i]))
	# num<-length(na.omit(Trait[,i]));
	# normt<-shapiro.test(na.omit(Trait[,i]))
	data<-rbind(data,data.frame(TraitName=TraitName[i],PhenotypeNum=num,Average=mean,Variance=var,W=normt$statistic,Pvalue=normt$p.value))
}
write.table(file=opt$outfile,data,quote=F,sep="\t",row.names=F);
