#!/usr/bin/env Rscript

times<-Sys.time()
library('getopt')
library(ggplot2)
library(grid)
library(gplots)

options(bitmapType='cairo')

spec = matrix(c(
	'GQ','a',0,'character',
	'dep','b',0,'character',
	'o','c',0,'character',
	'sub','d',0,'logical',
	'par','e',0,'character',
	'str','s',0,'character',
	'help' , 'f', 0, 'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript snp_qual.r --GQ  --dep  --o  
	
Usage:
	--GQ	Cumulative GQ
	--dep	Cumulative depth
	--o	    figure name
	--s	output variant type default SNP
	--help		usage
\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$GQ) ) { print_usage(spec) }
if ( is.null(opt$dep) ) { print_usage(spec) }
if ( is.null(opt$o) ) { print_usage(spec) }
if ( is.null(opt$str)){opt$srt="SNP"}

old_theme <- theme_update(
  axis.ticks=element_line(colour="black"),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  panel.background=element_blank(),
  panel.border=element_rect(fill="transparent", color="black"),
  axis.line=element_line(size=0.5),
  plot.title=element_text(hjust=0.5,size=10)
  
)

depth_in=read.table(opt$dep,header=T)
samples=unique(depth_in$sampleID);
for (i in samples){
	depth_sum=max(depth_in$num[depth_in$sampleID == i]);
	depth_in$frac[depth_in$sampleID == i]=depth_in$num[depth_in$sampleID == i]/depth_sum;
}
qual=read.table(opt$GQ,header=T)
samples=unique(qual$sampleID);
for (i in samples){ 
	qual_sum=max(qual$num[qual$sampleID == i]);
	qual$frac[qual$sampleID == i]=qual$num[qual$sampleID == i]/qual_sum;
}
pdf(paste(opt$o,".pdf",sep=""),width=12,height=6)
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,12)))
p<-ggplot(depth_in,aes(x=Depth,y=frac))+geom_line(aes(colour=sampleID,group=sampleID))+xlim(0,50)+theme(legend.position='none');
p<-p+xlab("Depth")+ylab(paste(opt$str," Fraction(%)",sep=""))+ggtitle(paste("Cumulative ",opt$str," depth 
distribution",sep=""))

print(p, vp = vplayout(1, 1:5.5))
p1<-ggplot(qual,aes(x=GQvalue,y=frac))+geom_line(aes(colour=sampleID,group=sampleID))
p1<-p1+xlab("Quality")+ylab(paste(opt$str," Fraction(%)",sep=""))+ggtitle(paste("Cumulative ",opt$str," depth 
distribution",sep=""))+guides(colour=guide_legend(title=NULL))

print(p1,vp = vplayout(1, 6:12))
dev.off();

png(paste(opt$o,".png",sep=""),width=1200,height=600)
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,1100)))
p<-ggplot(depth_in,aes(x=Depth,y=frac))+geom_line(aes(colour=sampleID,group=sampleID))+xlim(0,50)+theme(legend.position='none');
p<-p+xlab("Depth")+ylab(paste(opt$str," Fraction(%)",sep=""))+ggtitle(paste("Cumulative ",opt$str," depth 
distribution",sep=""))+guides(colour=guide_legend(title=NULL))

print(p, vp = vplayout(1, 1:500))
p1<-ggplot(qual,aes(x=GQvalue,y=frac))+geom_line(aes(colour=sampleID,group=sampleID))
p1<-p1+xlab("Quality")+ylab(paste(opt$str," Fraction(%)",sep=""))+ggtitle(paste("Cumulative ",opt$str," depth 
distribution",sep=""))+guides(colour=guide_legend(title=NULL))

print(p1,vp = vplayout(1, 500:1100))
dev.off();

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
