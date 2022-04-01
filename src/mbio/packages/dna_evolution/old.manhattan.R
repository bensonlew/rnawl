library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'fst','f',0,'character',
	'pi1','1',0,'character',
	'pi2','2',0,'character',
	'tajima1','3',0,'character',
	'tajima2','4',0,'character',
	'out','o',0,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--fst	the fst file
	--pi1	the pi1 file 
	--pi2	the pi2 file
	--tajima1 the tajima1 file
	--tajima2 the tajima2 file
	--out the output file
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$out)){ print_usage(spec) }
if ( is.null(opt$fst)){ print_usage(spec) }
if ( is.null(opt$pi1)){ print_usage(spec) }
if ( is.null(opt$pi2)){ print_usage(spec) }
if ( is.null(opt$tajima1)){ print_usage(spec) }
if ( is.null(opt$tajima2)){ print_usage(spec) }

library(qqman)
times<-Sys.time()
pi1<-read.table(opt$pi1,sep="\t",head=TRUE,stringsAsFactors=FALSE)
pi2<-read.table(opt$pi2,sep="\t",head=TRUE,stringsAsFactors=FALSE)
fst<-read.table(opt$fst,sep="\t",head=TRUE,stringsAsFactors=FALSE)
tajima1<-read.table(opt$tajima1,sep="\t",head=TRUE,stringsAsFactors=FALSE)
tajima2<-read.table(opt$tajima2,sep="\t",head=TRUE,stringsAsFactors=FALSE)
pi1$win<-paste(pi1$CHROM,pi1$BIN_START)
pi2$win<-paste(pi2$CHROM,pi2$BIN_START)
fst$win<-paste(fst$CHROM,fst$BIN_START)
tajima1$win<-paste(tajima1$CHROM,tajima1$BIN_START+1)
tajima2$win<-paste(tajima2$CHROM,tajima2$BIN_START+1)
win<-intersect(pi1$win,intersect(pi2$win,intersect(fst$win,intersect(tajima1$win,tajima2$win))))
data<-data.frame(chr=pi1$CHROM[pi1$win %in% win],stringsAsFactors=FALSE)
data$pos<-pi1$BIN_START[pi1$win %in% win]
data$pi1<-pi1$PI[pi1$win %in% win]
data$pi2<-pi2$PI[pi2$win %in% win]
data$tajima1<-tajima1$TajimaD[tajima1$win %in% win]
data$tajima2<-tajima2$TajimaD[tajima2$win %in% win]
data$fst<-fst$WEIGHTED_FST[fst$win %in% win]
data$tajima1[data$tajima1=="NaN"]=0
data$tajima2[data$tajima2=="NaN"]=0

write.table(file=paste(opt$out,".pi_tajimaD_fst.detail",sep=""),row.names=FALSE,data)
write.table(file=paste(opt$out,".pop1.pi_tajimaD_fst.select",sep=""),row.names=FALSE,subset(data,(data$pi1 < quantile(data$pi1,prob=0.05) & (data$tajima1 < quantile(data$tajima1,prob=0.05) | data$tajima1 > quantile(data$tajima1,prob=0.95)) & data$fst > quantile(data$fst,prob=0.95))))
write.table(file=paste(opt$out,".pop2.pi_tajimaD_fst.select",sep=""),row.names=FALSE,subset(data,(data$pi2 < quantile(data$pi2,prob=0.05) & (data$tajima2 < quantile(data$tajima2,prob=0.05) | data$tajima2 > quantile(data$tajima2,prob=0.95))& data$fst > quantile(data$fst,prob=0.95))))
#write.table(file=paste(opt$out,".diver.select",sep=""),row.names=FALSE,subset(data,(data$fst > quantile(data$fst,prob=0.95))))
#write.table(file=paste(opt$out,".fst.select",sep=""),row.names=FALSE,subset(data,(data$fst > quantile(data$fst,prob=0.95))))
chrlab=unique(data$chr)
for (i in 1:length(chrlab)){data$chr[data$chr==chrlab[i]]=i}
data$chr<-as.numeric(data$chr)
ylimit=max(abs(data$tajima1),abs(data$tajima2))
threpi1<-quantile(-log10(data$pi1),prob=0.95)
threpi2<-quantile(-log10(data$pi2),prob=0.95)
pdf(paste(opt$out,"pop1.manhattan.pdf",sep="."),height=12,width=16)
par(mfrow = c(3, 1))
manhattan(data,chr="chr",bp="pos",p="pi1",col=rainbow(4),chrlabs=chrlab,ylab="-log10(Pi)",logp=TRUE,suggestiveline=threpi1,genomewideline = F)
manhattan(data,chr="chr",bp="pos",p="tajima1",col=rainbow(4),chrlabs=chrlab,ylab="TajimaD",logp=FALSE,ylim=c(-1*ylimit,ylimit),suggestiveline=quantile(data$tajima1,prob=0.95),genomewideline = quantile(data$tajima1,prob=0.05))
manhattan(data,chr="chr",bp="pos",p="fst",col=rainbow(4),chrlabs=chrlab,ylab="Fst",logp=FALSE,suggestiveline=quantile(data$fst,prob=0.95),genomewideline = F)
dev.off()
pdf(paste(opt$out,"pop2.manhattan.pdf",sep="."),height=12,width=16)
par(mfrow = c(3, 1))
manhattan(data,chr="chr",bp="pos",p="pi2",col=rainbow(4),chrlabs=chrlab,ylab="-log10(Pi)",logp=TRUE,suggestiveline=threpi2,genomewideline = F)
manhattan(data,chr="chr",bp="pos",p="tajima2",col=rainbow(4),chrlabs=chrlab,ylab="TajimaD",logp=FALSE,,ylim=c(-1*ylimit,ylimit),suggestiveline=quantile(data$tajima2,prob=0.95),genomewideline = quantile(data$tajima2,prob=0.05))
manhattan(data,chr="chr",bp="pos",p="fst",col=rainbow(4),chrlabs=chrlab,ylab="Fst",ylim=c(0,1),logp=FALSE,suggestiveline=quantile(data$fst,prob=0.95),genomewideline = F)
dev.off()
png(paste(opt$out,"pop1.manhattan.png",sep="."),height=1200,width=1600)
par(mfrow = c(3, 1))
manhattan(data,chr="chr",bp="pos",p="pi1",col=rainbow(4),chrlabs=chrlab,ylab="-log10(Pi)",logp=TRUE,suggestiveline=threpi1,genomewideline = F)
manhattan(data,chr="chr",bp="pos",p="tajima1",col=rainbow(4),chrlabs=chrlab,ylab="TajimaD",logp=FALSE,ylim=c(-1*ylimit,ylimit),suggestiveline=quantile(data$tajima1,prob=0.95),genomewideline = quantile(data$tajima1,prob=0.05))
manhattan(data,chr="chr",bp="pos",p="fst",col=rainbow(4),chrlabs=chrlab,ylab="Fst",logp=FALSE,suggestiveline=quantile(data$fst,prob=0.95),genomewideline = F)
dev.off()
png(paste(opt$out,"pop2.manhattan.png",sep="."),height=1200,width=1600)
par(mfrow = c(3, 1))
manhattan(data,chr="chr",bp="pos",p="pi2",col=rainbow(4),chrlabs=chrlab,ylab="-log10(Pi)",logp=TRUE,suggestiveline=threpi2,genomewideline = F)
manhattan(data,chr="chr",bp="pos",p="tajima2",col=rainbow(4),chrlabs=chrlab,ylab="TajimaD",logp=FALSE,,ylim=c(-1*ylimit,ylimit),suggestiveline=quantile(data$tajima2,prob=0.95),genomewideline = quantile(data$tajima2,prob=0.05))
manhattan(data,chr="chr",bp="pos",p="fst",col=rainbow(4),chrlabs=chrlab,ylab="Fst",ylim=c(0,1),logp=FALSE,suggestiveline=quantile(data$fst,prob=0.95),genomewideline = F)
dev.off()


escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
