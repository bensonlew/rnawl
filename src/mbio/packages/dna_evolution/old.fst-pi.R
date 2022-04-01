library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'fst','f',1,'character',
	'pi1','1',1,'character',
	'pi2','2',1,'character',
	'out','o',1,'character',
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

library(ggplot2)
library(grid)
library("gridExtra")

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
write.table(file=paste(opt$out,"fst_pi.detail.select",sep="."),row.names=FALSE,subset(pi,pi$theta < quantile(pi$theta,probs=0.05) & pi$fst > quantile(pi$fst,probs=0.95) ));

empty<-ggplot()+theme(panel.background=element_blank())
hist_top<-ggplot()+theme_bw()
#hist_top<-hist_top+geom_histogram(aes(pi$theta[pi$theta >= quantile(pi$theta,probs=0.95)]),colour='green',fill='green',binwidth=0.01)
hist_top<-hist_top+geom_histogram(aes(pi$theta[pi$theta <= quantile(pi$theta,probs=0.05)]),colour='blue',fill='blue',binwidth=0.01)
hist_top<-hist_top+geom_histogram(aes(pi$theta[pi$theta > quantile(pi$theta,probs=0.05) & pi$theta < quantile(pi$theta,probs=0.95)]),colour='gray',fill='gray',binwidth=0.01)
hist_top<-hist_top+theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank())+theme(panel.background=element_blank())+ylab("Pi count")
hist_top<-hist_top+geom_vline(aes(xintercept=quantile(pi$theta,probs=0.05)),linetype=5,col="black")
hist_right<-ggplot()+theme_bw()
hist_right<-hist_right+geom_histogram(aes(pi$fst[pi$fst >= quantile(pi$fst,probs=0.95)]),colour='orange',fill='orange',binwidth=0.01)
hist_right<-hist_right+geom_histogram(aes(pi$fst[pi$fst <= quantile(pi$fst,probs=0.95)]),colour='gray',fill='gray',binwidth=0.01)
hist_right<-hist_right+theme(axis.title.y=element_blank(),axis.text.x = element_text(size=5))+theme(panel.background=element_blank())+coord_flip()+xlab("Fst distribution")
hist_right<-hist_right+geom_vline(aes(xintercept=quantile(pi$fst,probs=0.95)),linetype=5,col="black")+ylab("Fst Count")
scatter<-ggplot()+theme_bw()
scatter<-scatter+geom_point(aes(pi$theta[pi$theta > quantile(pi$theta,probs=0.05) & pi$theta < quantile(pi$theta,probs=0.95)],pi$fst[pi$theta > quantile(pi$theta,probs=0.05) & pi$theta < quantile(pi$theta,probs=0.95)]),colour='gray',fill='gray')
scatter<-scatter+geom_point(aes(pi$theta[pi$theta < quantile(pi$theta,probs=0.05) & pi$fst < quantile(pi$fst,probs=0.95)],pi$fst[pi$theta < quantile(pi$theta,probs=0.05) & pi$fst < quantile(pi$fst,probs=0.95)]),colour='gray',fill='gray')
scatter<-scatter+geom_point(aes(pi$theta[pi$theta < quantile(pi$theta,probs=0.05) & pi$fst > quantile(pi$fst,probs=0.95)],pi$fst[pi$theta < quantile(pi$theta,probs=0.05) & pi$fst > quantile(pi$fst,probs=0.95)]),colour='blue',fill='blue')
scatter<-scatter+geom_point(aes(pi$theta[pi$theta > quantile(pi$theta,probs=0.95) & pi$fst < quantile(pi$fst,probs=0.95)],pi$fst[pi$theta > quantile(pi$theta,probs=0.95) & pi$fst < quantile(pi$fst,probs=0.95)]),colour='gray',fill='gray')
scatter<-scatter+geom_point(aes(pi$theta[pi$theta > quantile(pi$theta,probs=0.95) & pi$fst > quantile(pi$fst,probs=0.95)],pi$fst[pi$theta > quantile(pi$theta,probs=0.95) & pi$fst > quantile(pi$fst,probs=0.95)]),colour='gray',fill='gray')
scatter<-scatter+theme(panel.background=element_blank())+xlab("Pi Distribution")+ylab("Fst Distribution")
scatter<-scatter+geom_hline(aes(yintercept=quantile(pi$fst,probs=0.95)),linetype=5,col="black")
#scatter<-scatter+geom_vline(aes(xintercept=quantile(pi$theta,probs=0.95)),linetype=5,col="black")
scatter<-scatter+geom_vline(aes(xintercept=quantile(pi$theta,probs=0.05)),linetype=5,col="black")



pdf(paste(opt$out,"fst_pi.pdf",sep="."));
grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4,1), heights=c(1,4))
dev.off()
png(paste(opt$out,"fst_pi.png",sep="."));
grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4,1), heights=c(1,4))
dev.off()

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)