library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'pi','p',1,'character',
	'tajima','t',1,'character',
	'out','o',1,'character',
	'thre','r',0,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4)
opt = getopt(spec)
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage:
	--pi	pi1 file
	--tajima     tajima type
	--out	out file
	--thre  threshold default 0.05
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
if ( is.null(opt$thre)){thre=0.05}else{thre=opt$thre}

thre=as.numeric(thre)
thre1=thre
thre2=1-thre
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
write.table(file=paste(opt$out,"pi_tajimaD.detail",sep="."),row.names=FALSE,pi);
data<-pi
write.table(file=paste(opt$out,"pi_tajimaD.detail.select",sep="."),row.names=FALSE,subset(data,data$pi < quantile(data$pi,prob=thre1) & (data$tajima < quantile(data$pi,prob=thre1) | data$tajima > quantile(data$pi,prob=thre2))));
# empty<-ggplot()+theme(panel.background=element_blank())
#
# diffpi<-(max(pi$pi)-min(pi$pi))/2000
# theta1<-ceiling((quantile(pi$pi,thre1)-quantile(pi$pi,0))/diffpi)
#
# g1<-ggplot()+theme_bw()
# g1<-g1+geom_histogram(aes(pi$pi,y=..density..),col=c(rep("green",theta1),rep("grey",2001-theta1)),binwidth=diffpi)
# g1<-g1+theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank())+theme(panel.background=element_blank())+ylab("Pi Frequency %")
# g1<-g1+geom_vline(aes(xintercept=quantile(pi$pi,probs=thre1)),linetype=5,col="black")
#
# difft<-(max(pi$tajima)-min(pi$tajima))/2000
# tajima1<-ceiling((quantile(pi$tajima,1)-quantile(pi$tajima,thre2))/difft)
# tajima2<-ceiling((quantile(pi$tajima,thre1)-quantile(pi$tajima,0))/difft)
#
#
# g2=ggplot()+theme_bw()
# g2=g2+geom_histogram(aes(pi$tajima,y=..density..),colour=c(rep('green',tajima2),rep('grey',2001-tajima1-tajima2),rep('blue',tajima1)),binwidth=difft)
# g2=g2+theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+theme(panel.background=element_blank())+coord_flip()+xlab("Tajima'D Count")
# g2=g2+geom_vline(aes(xintercept=quantile(pi$tajima,probs=thre2)),linetype=5,col="black")+ylab("Tajima'D Frequency %")
# g2=g2+geom_vline(aes(xintercept=quantile(pi$tajima,probs=thre1)),linetype=5,col="black")
#
#
#
# scatter<-ggplot()+theme_bw()
# scatter<-scatter+geom_point(aes(pi$pi,pi$tajima),colour='gray',fill='gray')
# scatter<-scatter+geom_point(aes(pi$pi[pi$pi <= quantile(pi$pi,probs=thre1) & pi$tajima < quantile(pi$tajima,probs=thre1)],pi$tajima[pi$pi <= quantile(pi$pi,probs=thre1) & pi$tajima < quantile(pi$tajima,probs=thre1)]),colour='gray',fill='gray')
# scatter<-scatter+geom_point(aes(pi$pi[pi$pi <= quantile(pi$pi,probs=thre1) & pi$tajima < quantile(pi$tajima,probs=thre1)],pi$tajima[pi$pi <= quantile(pi$pi,probs=thre1) & pi$tajima < quantile(pi$tajima,probs=thre1)]),colour='blue',fill='blue')
# scatter<-scatter+geom_point(aes(pi$pi[pi$pi <= quantile(pi$pi,probs=thre1) & pi$tajima > quantile(pi$tajima,probs=thre2)],pi$tajima[pi$pi <= quantile(pi$pi,probs=thre1) & pi$tajima > quantile(pi$tajima,probs=thre2)]),colour='green',fill='green')
# scatter<-scatter+theme(panel.background=element_blank())+xlab("Pi")+ylab("Tajima'D")
# scatter<-scatter+geom_hline(aes(yintercept=quantile(pi$tajima,probs=thre2)),linetype=5,col="black")
# scatter<-scatter+geom_hline(aes(yintercept=quantile(pi$tajima,probs=thre1)),linetype=5,col="black")
# scatter<-scatter+geom_vline(aes(xintercept=quantile(pi$pi,probs=thre1)),linetype=5,col="black")
#
#
# pdf(paste(opt$out,"pi_tajimaD.pdf",sep="."));
# # ggarrange(g1,empty,scatter,g2,widths=c(3,1),heights=c(1,3))
# dev.off()
# png(paste(opt$out,"pi_tajimaD.png",sep="."));
# # ggarrange(g1,empty,scatter,g2,widths=c(3,1),heights=c(1,3))
# dev.off()

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
