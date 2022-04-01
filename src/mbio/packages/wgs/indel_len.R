#!/usr/bin/env Rscript

times<-Sys.time()
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'i','a',0,'character',
	'od','b',0,'character',
	'help','c',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
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
df<-read.table(opt$i,header=T)
sample<-unique(df[,1])

for (i in 1:length(sample)){
	options(scipen=200)
	pdf(file=paste(opt$o,"/",sample[i],".indel.pdf",sep=""))
	a=df[df[,1]%in%sample[i],2]
	xlim_max_num=10
	xlim_min_num=-10
	bins=seq(xlim_min_num,xlim_max_num,by=1)
	 hist(a[a>=-10 & a<=10],breaks=bins,col="blue",xlab="Indel Lenth (bp)",ylab="number",freq=TRUE,main="Indel size distribution",xaxt="n",xlim=c(xlim_min_num,xlim_max_num))
	 axis= axis(1,at=seq(xlim_min_num,xlim_max_num,5))
	box()
	dev.off()
	png(file=paste(opt$o,"/",sample[i],".indel.png",sep=""))
	a=df[df[,1]%in%sample[i],2]
	hist(a[a>=-10 & a<=10],breaks=bins,col="blue",xlab="Indel Lenth (bp)",ylab="number",freq=TRUE,main="Indel size distribution",xaxt="n",xlim=c(xlim_min_num,xlim_max_num))
	axis= axis(1,at=seq(xlim_min_num,xlim_max_num,5))
	box()
	dev.off()

}
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)





#TSHKO.insertsize.pdf
