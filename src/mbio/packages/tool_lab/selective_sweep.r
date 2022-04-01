library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'fst','f',1,'character',
	'pi1','1',1,'character',
	'pi2','2',1,'character',
	'out','o',1,'character',
	'thre','t',0,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4)
opt = getopt(spec)
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage:
	--fst	fst file
	--pi1	pi1 file
	--pi2	pi2 file
	--thre  threshold default 0.05
	--out	out file
	--help		usage
\n")
	q(status=1);
}
times<-Sys.time()
library('qtl');
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$fst) ) { print_usage(spec) }
if ( is.null(opt$pi1) ) { print_usage(spec) }
if ( is.null(opt$pi2) ) { print_usage(spec) }
if ( is.null(opt$out) ) { opt$out="./";}
if ( is.null(opt$thre)){thre=0.05}else{thre=opt$thre}

thre=as.numeric(thre)
thre1=thre
thre2=1-thre

library(ggplot2)
library(grid)
library("gridExtra")
library(egg)
pi1<-read.table(opt$pi1,sep="\t",head=TRUE)
pi2<-read.table(opt$pi2,sep="\t",head=TRUE)
fst<-read.table(opt$fst,sep="\t",head=TRUE)

pi1$win<-paste(pi1$CHROM,pi1$BIN_START,pi1$BIN_END)
pi2$win<-paste(pi2$CHROM,pi2$BIN_START,pi2$BIN_END)
fst$win<-paste(fst$CHROM,fst$BIN_START,fst$BIN_END)
win<-intersect(intersect(pi1$win,pi2$win),fst$win)
pi<-data.frame(chr=pi1$CHROM[pi1$win %in% win],pos1=pi1$BIN_START[pi1$win %in% win],pos2=pi1$BIN_END[pi1$win %in% win],pi1=pi1$PI[pi1$win %in% win],pi2=pi2$PI[pi2$win %in% win],fst=fst$WEIGHTED_FST[fst$win %in% win])
pi$theta=pi$pi1/pi$pi2
pi$fst[pi$fst < 0]=0
write.table(file=paste(opt$out,"fst_pi.detail",sep="."),row.names=FALSE,pi);
write.table(file=paste(opt$out,"fst_pi.detail-pi1.select",sep="."),row.names=FALSE,subset(pi,pi$theta < quantile(pi$theta,probs=thre1) & pi$fst > quantile(pi$fst,probs=thre2) ))
write.table(file=paste(opt$out,"fst_pi.detail-pi2.select",sep="."),row.names=FALSE,subset(pi,pi$theta > quantile(pi$theta,probs=thre2) & pi$fst > quantile(pi$fst,probs=thre2) ))

empty<-ggplot()+theme(panel.background=element_blank())

diff<-(max(pi$theta)-min(pi$theta))/2000

theta1<-ceiling((quantile(pi$theta,thre1)-quantile(pi$theta,0))/diff)
print(quantile(pi$theta,0))
theta2<-ceiling((quantile(pi$theta,1)-quantile(pi$theta,thre2))/diff)
g1<-ggplot()+theme_bw()
g1<-g1+geom_histogram(aes(pi$theta,y=..density..),col=c(rep('blue',theta1),rep('grey',2001-theta1-theta2),rep('green',theta2)),binwidth=diff)
g1<-g1+theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank())+theme(panel.background=element_blank())+ylab("Pi Ratio Frequency %")
g1<-g1+geom_vline(aes(xintercept=quantile(pi$theta,probs=thre2)),linetype=5,col="black")
g1<-g1+geom_vline(aes(xintercept=quantile(pi$theta,probs=thre1)),linetype=5,col="black")
scatter<-ggplot()+theme_bw()
scatter<-scatter+geom_point(aes(pi$theta[pi$theta > quantile(pi$theta,probs=thre1) & pi$theta < quantile(pi$theta,probs=thre2)],pi$fst[pi$theta > quantile(pi$theta,probs=thre1) & pi$theta < quantile(pi$theta,probs=thre2)]),colour='gray',fill='gray')
scatter<-scatter+geom_point(aes(pi$theta[pi$theta < quantile(pi$theta,probs=thre1) & pi$fst < quantile(pi$fst,probs=thre2)],pi$fst[pi$theta < quantile(pi$theta,probs=thre1) & pi$fst < quantile(pi$fst,probs=thre2)]),colour='gray',fill='gray')
scatter<-scatter+geom_point(aes(pi$theta[pi$theta < quantile(pi$theta,probs=thre1) & pi$fst > quantile(pi$fst,probs=thre2)],pi$fst[pi$theta < quantile(pi$theta,probs=thre1) & pi$fst > quantile(pi$fst,probs=thre2)]),colour='blue',fill='blue')
scatter<-scatter+geom_point(aes(pi$theta[pi$theta > quantile(pi$theta,probs=thre2) & pi$fst < quantile(pi$fst,probs=thre2)],pi$fst[pi$theta > quantile(pi$theta,probs=thre2) & pi$fst < quantile(pi$fst,probs=thre2)]),colour='gray',fill='gray')
scatter<-scatter+geom_point(aes(pi$theta[pi$theta > quantile(pi$theta,probs=thre2) & pi$fst > quantile(pi$fst,probs=thre2)],pi$fst[pi$theta > quantile(pi$theta,probs=thre2) & pi$fst > quantile(pi$fst,probs=thre2)]),colour='green',fill='green')
scatter<-scatter+theme(panel.background=element_blank())+xlab("Pi Ratio")+ylab("Fst")
scatter<-scatter+geom_hline(aes(yintercept=quantile(pi$fst,probs=thre2)),linetype=5,col="black")
scatter<-scatter+geom_vline(aes(xintercept=quantile(pi$theta,probs=thre2)),linetype=5,col="black")
scatter<-scatter+geom_vline(aes(xintercept=quantile(pi$theta,probs=thre1)),linetype=5,col="black")
difffst<-(max(pi$fst)-min(pi$fst))/2000
fst1<-ceiling((quantile(pi$fst,1)-quantile(pi$fst,thre2))/difffst)
g3<-ggplot()+theme_bw()
g3<-g3+geom_histogram(aes(pi$fst,y=..density..),colour=c(rep('grey',2001-fst1),rep('orange',fst1)),binwidth=difffst)+ylim(0,8)
g3<-g3+theme(axis.title.y=element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())+theme(panel.background=element_blank())+coord_flip()+xlab("Fst distribution")
g3<-g3+geom_vline(aes(xintercept=quantile(pi$fst,probs=thre2)),linetype=5,col="black")+ylab("Fst Frequency %")

pdf(paste(opt$out,"fst_pi.pdf",sep="."), onefile=FALSE);
ggarrange(g1,empty,scatter,g3,widths=c(3,1),heights=c(1,3))
dev.off()
png(paste(opt$out,"fst_pi.png",sep="."));
ggarrange(g1,empty,scatter,g3,widths=c(3,1),heights=c(1,3))
dev.off()




escaptime=Sys.time()-times;
print("Done!")
print(escaptime)


