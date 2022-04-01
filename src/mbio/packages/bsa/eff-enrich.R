#!/usr/bin/env Rscript
library('getopt');
library(ggplot2)
library(grid)

options(bitmapType='cairo')
spec = matrix(c(
	'input','i',0,'character',
	'output','o',0,'character',
	'top','t',0,'character',
	'eggnog','e',0,'character',
	'help','m',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage:
	--input	the input  file
	--output	the out file
	--top the top file
	--eggnog the eggnog for draw
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$input)) { print_usage(spec)}
if ( is.null(opt$output)){ print_usage(spec) }
times<-Sys.time()

data<-read.table(opt$input,sep="\t",comment.char="^",head=TRUE)
names(data)<-c("id","description","k","M","n","N");
data<-na.omit(data)
pvalue<-phyper(data$k,data$M,data$N-data$M,data$n,lower.tail=FALSE);
qvalue<-p.adjust(pvalue,method="fdr")
outd<-data.frame(id=data$id,des=data$description,eff=data$k,total=data$M,pvalue=1-pvalue,qvalue=1-qvalue)
outd<-outd[order(outd$eff,decreasing=T),]
write.table(file=paste(opt$output,"detail",sep="."),outd,row.name=FALSE,sep="\t");
# draw<-data.frame(id=outd$id,des=outd$des,eff=outd$eff)
# if(!is.null(opt$top)){
# 	draw<-data.frame(id=outd$id[1:20],des=outd$des[1:20],eff=outd$eff[1:20])
# }
# draw$labels=paste(draw$id,draw$des,sep=":")
# if(is.null(opt$eggnog)){
# 	p<-ggplot(draw,aes(x=labels,y=eff))+ theme_bw()
# 	p<-p+geom_bar(stat="identity",aes(fill=labels),width=0.5)
# 	p<-p+theme(plot.title=element_text(hjust=0.5),legend.position="NONE")
# 	p<-p+labs(x="Function", title="Most Effdata of Gene Function",y="Eff Number")+coord_flip()
# }else{
# 	p<-ggplot(draw,aes(x=id,y=eff))+ theme_bw()
# 	p<-p+geom_bar(stat="identity",aes(fill=labels),width=0.5)+scale_fill_discrete(name="")
# 	p<-p+theme(plot.title=element_text(hjust=0.5),legend.position="right",legend.key.size=unit(0.8,"cm"))
# 	p<-p+labs(x="Function", title="Most Effdata of Gene Function",y="Eff Number")+ guides(fill = guide_legend(ncol = 1))
# }
# ggsave(filename=paste(opt$output,"pdf",sep="."),p,height=9,width=16,device="pdf")
# ggsave(filename=paste(opt$output,"png",sep="."),p,height=9,width=16,device="png")


escaptime=Sys.time()-times;

print("Done!")
print(escaptime)
