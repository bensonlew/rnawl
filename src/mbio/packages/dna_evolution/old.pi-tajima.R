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
pi<-read.table(opt$pi,sep="\t",head=TRUE)
tajima<-read.table(opt$tajima,sep="\t",head=TRUE)
pi$win<-paste(pi$CHROM,pi$BIN_START)
tajima$win<-paste(tajima$CHROM,tajima$BIN_START+1)
win<-intersect(tajima$win,pi$win)

pi<-data.frame(chr=pi$CHROM[pi$win %in% win],pos1=pi$BIN_START[pi$win %in% win],pos2=pi$BIN_END[pi$win %in% win],pi=pi$PI[pi$win %in% win],tajima=tajima$TajimaD[tajima$win %in% win])
pi$tajima[pi$tajima == "NaN"]=0
write.table(file=paste(opt$out,"pi_tajimaD.detail",sep="."),row.names=FALSE,pi);
data<-pi
write.table(file=paste(opt$out,"pi_tajimaD.detail.select",sep="."),row.names=FALSE,subset(data,data$pi < quantile(data$pi,prob=0.05) & (data$tajima < quantile(data$pi,prob=0.05) | data$tajima > quantile(data$pi,prob=0.95))));
empty<-ggplot()+theme(panel.background=element_blank())
maxpi=max(pi$pi)
minpi=min(pi$pi)
hist_top<-ggplot()+theme_bw()
hist_top<-hist_top+geom_histogram(aes(pi$pi),colour='gray',fill='gray',binwidth=(maxpi-minpi)/1000)
hist_top<-hist_top+geom_histogram(aes(pi$pi[pi$pi <= quantile(pi$pi,probs=0.05)]),colour='blue',fill='blue',binwidth=(maxpi-minpi)/1000)
hist_top<-hist_top+theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank())+theme(panel.background=element_blank())+ylab("Pi Count")
hist_top<-hist_top+geom_vline(aes(xintercept=quantile(pi$pi,probs=0.05)),linetype=5,col="black")
hist_right=ggplot()+theme_bw()
hist_right=hist_right+geom_histogram(aes(pi$tajima[pi$tajima >= quantile(pi$tajima,probs=0.95)]),colour='green',fill='green',binwidth=0.01)
hist_right=hist_right+geom_histogram(aes(pi$tajima[pi$tajima <= quantile(pi$tajima,probs=0.05)]),colour='blue',fill='blue',binwidth=0.01)
hist_right=hist_right+geom_histogram(aes(pi$tajima[pi$tajima >= quantile(pi$tajima,probs=0.05) & pi$tajima <= quantile(pi$tajima,probs=0.95) ]),colour='gray',fill='gray',binwidth=0.01)
hist_right=hist_right+theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size=5))+theme(panel.background=element_blank())+coord_flip()+xlab("Tajima'D Distribution")
hist_right=hist_right+geom_vline(aes(xintercept=quantile(pi$tajima,probs=0.95)),linetype=5,col="black")+ylab("Tajima'D Count")
hist_right=hist_right+geom_vline(aes(xintercept=quantile(pi$tajima,probs=0.05)),linetype=5,col="black")
scatter<-ggplot()+theme_bw()
scatter<-scatter+geom_point(aes(pi$pi,pi$tajima),colour='gray',fill='gray')
scatter<-scatter+geom_point(aes(pi$pi[pi$pi <= quantile(pi$pi,probs=0.05) & pi$tajima < quantile(pi$tajima,probs=0.05)],pi$tajima[pi$pi <= quantile(pi$pi,probs=0.05) & pi$tajima < quantile(pi$tajima,probs=0.05)]),colour='gray',fill='gray')
scatter<-scatter+geom_point(aes(pi$pi[pi$pi <= quantile(pi$pi,probs=0.05) & pi$tajima < quantile(pi$tajima,probs=0.05)],pi$tajima[pi$pi <= quantile(pi$pi,probs=0.05) & pi$tajima < quantile(pi$tajima,probs=0.05)]),colour='blue',fill='blue')
scatter<-scatter+geom_point(aes(pi$pi[pi$pi <= quantile(pi$pi,probs=0.05) & pi$tajima > quantile(pi$tajima,probs=0.95)],pi$tajima[pi$pi <= quantile(pi$pi,probs=0.05) & pi$tajima > quantile(pi$tajima,probs=0.95)]),colour='green',fill='green')
scatter<-scatter+theme(panel.background=element_blank())+xlab("Pi Distribution")+ylab("Tajima'D Distribution")
scatter<-scatter+geom_hline(aes(yintercept=quantile(pi$tajima,probs=0.95)),linetype=5,col="black")
scatter<-scatter+geom_hline(aes(yintercept=quantile(pi$tajima,probs=0.05)),linetype=5,col="black")
scatter<-scatter+geom_vline(aes(xintercept=quantile(pi$pi,probs=0.05)),linetype=5,col="black")


pdf(paste(opt$out,"pi_tajimaD.pdf",sep="."));
grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4,1), heights=c(1,4))
dev.off()
png(paste(opt$out,"pi_tajimaD.png",sep="."));
grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4,1), heights=c(1,4))
dev.off()

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)