times<-Sys.time()
library('getopt');
options(bitmapType='cairo')
options(scipen = 500)
spec = matrix(c(
	'infile','i',0,'character',
	'outfile','o',0,'character',
	'col','c',0,'character',
	'win','w',0,'character',
	'step','s',0,'character',
	'method','m',0,'character',
	'abs','a',0,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage:
	--infile	the input hapmap file
	--outfile	the trait file
	--col	the col of chr pos index
	--win	the window size
	--step	the step size
	--method	default bp or num
	--abs	default 0 no abs adjust
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$infile)) { print_usage(spec)}
if ( is.null(opt$outfile)){ print_usage(spec) }
if ( is.null(opt$win)){opt$win=2000000;}
if ( is.null(opt$step)){opt$step=opt$win/40}
if ( is.null(opt$method)){opt$method="bp"}
if (is.null(opt$abs))(opt$abs="0");

data<-read.table(opt$infile,head=TRUE,comment.char = "^");
opt$col=length(unlist(strsplit(opt$col,",")))
if(opt$col ==3){
	chr<-data$X.chr
	pos<-data$pos
	index<-data$index
}else{
	chr<-data$X.chr
	pos<-data$pos
	index1<-data$INDEX1
	index2<-data$INDEX2
	delta<-data$DELTA
	index=delta
}
ewin<-function(pos,index,p1,p2){
	a<-index[pos > p1 & pos < p2 & index >= thres]
	return(length(a))
}
twin<-function(pos,p1,p2){
	a<-pos[pos > p1 & pos < p2]
	return(length(a))
}
mwin<-function(pos,index,p1,p2){
	a<-index[pos > p1 & pos < p2]
	return(mean(a))
}
print("sliding window")
chrname=unique(chr)
slid<-NULL
for (i in 1:(length(chrname))){
	print(chrname[i]);
	chrpos=pos[which(chr==chrname[i])]
	backpos=chrpos;
	if (opt$method=="num"){chrpos=c(1:length(chrpos))}
	if(opt$col ==3){
		chrindex=index[which(chr==chrname[i])]
	}else{
		chrindex1=index1[which(chr==chrname[i])]
		chrindex2=index2[which(chr==chrname[i])]
		chrdelta=delta[which(chr==chrname[i])]
	}
	chrlen=max(chrpos);
	win=ceiling(chrlen/as.numeric(opt$step));
	ss=c(1:win)
	pos1=ss*as.numeric(opt$step)-as.numeric(opt$win)/2;
	pos2=ss*as.numeric(opt$step)+as.numeric(opt$win)/2;
	if (opt$method=="num"){
		pos1[pos1 <= 1]=1;
		pos2[pos2 >= chrlen]=chrlen;
		pos1=backpos[pos1];
		pos2=backpos[pos2];
	}else{
		pos1[pos1 < 0]=0;
		pos2[pos2 > chrlen]=chrlen;
	}
	x=data.frame(pos1,pos2);
	if(opt$col ==3){
		wmean=apply(x,MARGIN=1,function(x,y,z,a) mwin(backpos,chrindex,x[1],x[2]));
		total=apply(x,MARGIN=1,function(x,y,z) twin(backpos,x[1],x[2]));
		slid<-rbind(slid,data.frame(chr=chrname[i],pos1=pos1,pos2=pos2,index=wmean,twin=total))
	}else{
		wmean1=apply(x,MARGIN=1,function(x,y,z,a) mwin(backpos,chrindex1,x[1],x[2]));
		wmean2=apply(x,MARGIN=1,function(x,y,z,a) mwin(backpos,chrindex2,x[1],x[2]));
		indexs=wmean1-wmean2;
		if(opt$abs == "1"){indexs=abs(indexs)}
		total=apply(x,MARGIN=1,function(x,y,z) twin(backpos,x[1],x[2]));
		slid<-rbind(slid,data.frame(chr=chrname[i],pos1=pos1,pos2=pos2,index1=wmean1,index2=wmean2,delta=indexs,twin=total))
	}
}
write.table(file=paste(opt$outfile,".slid.origin.result",sep=""),slid,row.name=FALSE,sep="\t")

total<-length(chr);
N=total
if(opt$col ==3){
	slid$index[slid$twin < 10]=mean(slid$index[slid$twin >= 10])
	thres<-quantile(slid$index[slid$twin > 10],probs=0.999,na.rm=TRUE)
	M=length(slid$index[slid$index > thres & slid$twin > 10])
}else{
	slid$index1[slid$twin < 10]=mean(slid$index1[slid$twin >= 10])
	slid$index2[slid$twin < 10]=mean(slid$index2[slid$twin >= 10])
	slid$delta[slid$twin < 10]=mean(slid$delta[slid$twin >= 10])
	thres<-quantile(slid$delta,probs=0.999,na.rm=TRUE)
	M=length(slid$delta[slid$delta > thres])
}
write.table(file=paste(opt$outfile,".slid.result",sep=""),slid,row.name=FALSE,sep="\t")

info<-NULL
for (i in 1:(length(chrname))){
	chrpos=pos[which(chr==chrname[i])]
	backpos=chrpos;
	if (opt$method=="num"){chrpos=c(1:length(chrpos))}
	#backpos=chrpos;
	chrindex=index[which(chr==chrname[i])]
	chrlen=max(chrpos);
	win=ceiling(chrlen/as.numeric(opt$step));
	ss=c(1:win)
	pos1=ss*as.numeric(opt$step)-as.numeric(opt$win)/2;
	pos2=ss*as.numeric(opt$step)+as.numeric(opt$win)/2;
	if (opt$method=="num"){
		pos1[pos1 <=1]=1;
		pos2[pos2 > chrlen]=chrlen;
		pos1=backpos[pos1]
		pos2=backpos[pos2]
	}else{
		pos1[pos1 < 0]=0;
		pos2[pos2 > chrlen]=chrlen;
	}
	x=data.frame(pos1,pos2);
	k=apply(x,MARGIN=1,function(x,y,z,a) ewin(backpos,chrindex,x[1],x[2]));
	n=apply(x,MARGIN=1,function(x,y,z) twin(backpos,x[1],x[2]));
	pvalue=phyper(k,M,N-M,n,lower.tail=FALSE);
	info<-rbind(info,data.frame(chr=chrname[i],pos1=pos1,pos2=pos2,k=k,M=M,N=N,n=n,Pvalue=pvalue))
}
fdr=p.adjust(info$Pvalue,method="bonferroni")

print(length(unique(info$chr)));

if(opt$col ==3){
	df=data.frame(chr=info$chr,pos1=info$pos1,pos2=info$pos2,index=slid$index,threhold=rep(thres,length(info$pos1)),total=slid$twin,peak=info$k,pvalue=info$Pvalue,fdr=fdr,stringsAsFactors=FALSE);
}else{
	df=data.frame(chr=info$chr,pos1=info$pos1,pos2=info$pos2,index1=slid$index1,index2=slid$index2,delta=slid$delta,threhold=rep(thres,length(info$pos1)),total=slid$twin,peak=info$k,pvalue=info$Pvalue,fdr=fdr,stringsAsFactors=FALSE);
}
df<-na.omit(df)
write.table(file=paste(opt$outfile,".result",sep=""),df,row.name=FALSE,sep="\t")

if(opt$col == 3){
	df$index[df$total < 10]=mean(df$index[df$total >= 10])
	write.table(file=paste(opt$outfile,".threshold.select",sep=""),subset(df,df$index > df$threhold),row.names=FALSE,sep="\t")
	write.table(file=paste(opt$outfile,".fdr.select",sep=""),subset(df,df$fdr < 0.01/N),row.names=FALSE,sep="\t")
}else{
	df$index1[df$total < 10]=mean(df$index1[df$total >= 10])
	df$index2[df$total < 10]=mean(df$index2[df$total >= 10])
	df$delta[df$total < 10]=mean(df$delta[df$total >= 10])
	write.table(file=paste(opt$outfile,".threshold.select",sep=""),subset(df,df$delta > df$threhold),row.names=FALSE,sep="\t")
	write.table(file=paste(opt$outfile,".fdr.select",sep=""),subset(df,df$delta < 0.01/N),row.names=FALSE,sep="\t")
}

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()
