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
dfpos<-data.frame(CHR=pos$X.chr,BP=pos$pos,Gprime=pos$Gprime,G=pos$G,Gfdr=pos$Gpval)
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
		p2 <- ggplot(dpos) +
		    geom_point(aes(x=BPcum, y=G,color=as.factor(order))) +geom_line(mapping = aes(x=BPcum,y=Gprime),color="black")+
		    scale_color_manual(values = rep(c("grey", "skyblue"), length(levels(dpos$CHR))))+
		    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
		    #geom_hline(yintercept = 0.9995,linetype="dashed",col="red")+
		    theme_bw() +xlab("chromosome")+ylab("G value")+
		    theme( 
		      legend.position="none",
		      panel.border = element_blank(),
		      panel.grid.major.x = element_blank(),
		      panel.grid.minor.x = element_blank()
		    )
		ggsave(file=paste(opt$output,"G.pdf",sep="."),p2,device="pdf",dpi=300,height=9,width=16)
		ggsave(file=paste(opt$output,"G.png",sep="."),p2,device="png",dpi=300,height=9,width=16)
		chrdata=filter(dpos,grepl(x=CHR,"chr"))
		chrdata$Gfdr=-1*log10(chrdata$Gfdr)
		p2 <- ggplot(chrdata) +
		    geom_point(aes(x=BP, y=Gfdr,color=as.factor(order))) +geom_hline(yintercept = -1*log10(0.05),linetype="dashed",col="red") +
		    scale_color_manual(values = rep(c("grey", "skyblue"), length(levels(dpos$CHR))))+
		    theme_bw() +xlab("chromosome")+ylab("-log10(p value)")+
		    theme( 
		      legend.position="none",
		      panel.border = element_blank(),
		      panel.grid.major.x = element_blank(),
		      panel.grid.minor.x = element_blank()
		    )
		    if(!is.null(opt$threshold)){
			p1=p1+geom_hline(yintercept = as.numeric(opt$threshold),linetype="dashed",col="red")
		    }
		ggsave(file=paste(opt$output,"pvalue.pdf",sep="."),p2,device="pdf",dpi=300,height=9,width=16)
		ggsave(file=paste(opt$output,"pvalue.png",sep="."),p2,device="png",dpi=300,height=9,width=16)


escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
