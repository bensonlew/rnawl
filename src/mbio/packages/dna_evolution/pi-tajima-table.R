library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'pi','p',1,'character',
	'tajima','t',1,'character',
	'out','o',1,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4)
opt = getopt(spec)
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example:

Usage:
	--pi	pi1 file
	--tajima	tajima type
	--out	out file
	--help		usage
\n")
	q(status=1);
}
times<-Sys.time()
library('qtl');
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$tajima) ) { print_usage(spec) }
if ( is.null(opt$pi) ) { print_usage(spec) }
if ( is.null(opt$out) ) { opt$out="./";}

library(ggplot2)
library(grid)
library("gridExtra")
# library(egg)
pi<-read.table(opt$pi,sep="\t",head=TRUE)
tajima<-read.table(opt$tajima,sep="\t",head=TRUE)
pi$win<-paste(pi$CHROM,pi$BIN_START)
tajima$win<-paste(tajima$CHROM,tajima$BIN_START+1)
win<-intersect(tajima$win,pi$win)

pi<-data.frame(chr=pi$CHROM[pi$win %in% win],pos1=pi$BIN_START[pi$win %in% win],pos2=pi$BIN_END[pi$win %in% win],pi=pi$PI[pi$win %in% win],tajima=tajima$TajimaD[tajima$win %in% win])
pi$tajima[pi$tajima == "NaN"]=0
tajima.table<-table(cut(pi$tajima,2000))
tajima<-as.vector(tajima.table)
tajima<tajima/sum(tajima)
tajima.rownames<-rownames(tajima.table)
df.tajima<-data.frame(x=tajima.rownames,y=tajima);
pi.table<-table(cut(pi$pi,2000))
pi<-as.vector(pi.table)
pi<pi/sum(pi)
pi.rownames<-rownames(pi.table)
df.pi<-data.frame(x=pi.rownames,y=pi);
write.table(file=paste(opt$out,"tajima.hist",sep="."),row.names=FALSE,df.tajima)
write.table(file=paste(opt$out,"pi.hist",sep="."),row.names=FALSE,df.pi)
write.table(file=paste(opt$out,"pi-tajima.detail",sep="."),row.names=FALSE,pi);




escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
