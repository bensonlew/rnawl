library('getopt');
options(bitmapType='cairo')
options(scipen = 200)
spec = matrix(c(
	'result','r',0,'character',
	'threshold','t',0,'character',
	'method','m',0,'character',
	'output','o',0,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--result	 the input snp file
	--method	method [G index ED]
	--thresold	the threshold default 0.9995
	--mutmap	the mutmap 
	--output	the out file 
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$result)) { print_usage(spec)}
if ( is.null(opt$output)){ opt$output="manhattan"}
if ( is.null(opt$method)){ opt$method="index" }
times<-Sys.time()
library(ggplot2)
library(dplyr)

pos<-read.table(opt$result,head=TRUE,comment.char="^")
pos$pos=(pos$pos1+pos$pos2)
dfpos<-data.frame(CHR=pos$X.chr,BP=pos$pos,Delta=pos$delta,Slide=pos$slidingD,CI=pos$CI)
#Gprime=pos$Gprime,G=pos$G,ED=pos$ED,EDprime=pos$EDprime,Gfdr=pos$Gfdr,EDfdr=pos$EDfdr)
	lev<-NULL
	lev$CHR<-levels(dfpos$CHR)
	lev$order<-gsub("chr","",lev$CHR)
	lev$order<-gsub("sca","1000",lev$order)
	lev$order=as.numeric(lev$order)
	dfpos=merge(dfpos,lev,by="CHR")
	dfpos=arrange(dfpos,order,BP)
	dpos <- dfpos %>% group_by(order) %>% summarise(chr_len=max(BP)) %>% mutate(tot=cumsum(chr_len)-chr_len) %>% select(-chr_len) %>%
	  left_join(dfpos, ., by=c("order"="order")) %>%
	  arrange(order, BP) %>%
	  mutate( BPcum=BP+tot)
	axisdf <- dpos %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
		p1 <- ggplot(dpos) +
		    geom_point(aes(x=BPcum, y=Delta,color=as.factor(order))) +geom_line(mapping = aes(x=BPcum,y=Slide),color="black")+
		    geom_line(mapping = aes(x=BPcum,y=CI),color="red")+
		    scale_color_manual(values = rep(c("grey", "skyblue"), length(levels(dpos$CHR))))+
		    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
		    theme_bw() +xlab("chromosome")+ylab("delta-index")+
		    theme( 
		      legend.position="none",
		      panel.border = element_blank(),
		      panel.grid.major.x = element_blank(),
		      panel.grid.minor.x = element_blank()
		    )
		    if(!is.null(opt$threshold)){
			p1=p1+geom_hline(yintercept = as.numeric(opt$threshold),linetype="dashed",col="red")
		    }
		ggsave(file=paste(opt$out,"index.pdf",sep="."),p1,device="pdf",dpi=300,height=9,width=16)
		ggsave(file=paste(opt$out,"index.png",sep="."),p1,device="png",dpi=300,height=9,width=16)
		chrdata=filter(dpos,grepl(x=CHR,"chr"))
		p1 <- ggplot(chrdata) +
		    geom_point(aes(x=BP, y=Delta,color=as.factor(order))) +geom_line(mapping = aes(x=BP,y=Slide),color="black")+
		    geom_line(mapping = aes(x=BP,y=CI),color="red")+
		    scale_color_manual(values = rep(c("grey", "skyblue"), length(levels(dpos$CHR))))+
		    theme_bw() +xlab("bp")+ylab("delta-index")+facet_wrap(~CHR,ncol=5)+
		    theme( 
		      legend.position="none",
		      panel.border = element_blank(),
		      panel.grid.major.x = element_blank(),
		      panel.grid.minor.x = element_blank()
		    )
		    if(!is.null(opt$threshold)){
			p1=p1+geom_hline(yintercept = as.numeric(opt$threshold),linetype="dashed",col="red")
		    }
		ggsave(file=paste(opt$out,"chr.index.pdf",sep="."),p1,device="pdf",dpi=300,height=9,width=16)
		ggsave(file=paste(opt$out,"chr.index.png",sep="."),p1,device="png",dpi=300,height=9,width=16)


escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
