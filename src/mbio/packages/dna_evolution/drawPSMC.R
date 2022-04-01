#!/usr/bin/env Rscript
library('getopt');
options(bitmapType='cairo')
options(scipen=200)
spec = matrix(c(
	'infile','i',0,'character',
	'outdir','o',0,'character',
	'year','y',0,'character',
	'help','m',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	q(status=1);
}
if ( !is.null(opt$help)) { print_usage(spec) }
if ( is.null(opt$infile)) { print_usage(spec)}
if ( is.null(opt$outdir)){ print_usage(spec) }
if(is.null(opt$year)){opt$year=1}
setwd(opt$outdir);
psmc.result<-function(file,i.iteration=25,mu=3.7e-8,s=100,g=as.numeric(opt$year))
{
	X<-scan(file=file,what="",sep="\n",quiet=TRUE)
	
	START<-grep("^RD",X)
	END<-grep("^//",X)
	
	X<-X[START[i.iteration+1]:END[i.iteration+1]]
	
	TR<-grep("^TR",X,value=TRUE)
	RS<-grep("^RS",X,value=TRUE)
	
	write(TR,"temp.psmc.result")
	theta0<-as.numeric(read.table("temp.psmc.result")[1,2])
	N0<-theta0/4/mu/s
	
	write(RS,"temp.psmc.result")
	a<-read.table("temp.psmc.result")
	Generation<-as.numeric(2*N0*a[,3])
	Ne<-as.numeric(N0*a[,4])
	
	file.remove("temp.psmc.result")
	
	n.points<-length(Ne)
	YearsAgo<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
		Generation[n.points])*g
	Ne<-c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
		Ne[n.points])
	
	data.frame(YearsAgo,Ne)
}
plotPsmc.popwise<-function(keywords, labels,
	command=NA, N0=25000,
	save.as="png", height=7, width=12,
	i.iteration=25, mu=3.7e-8, s=100, g=2,
	ylim=c(0,100000), xlim=c(200,500000),
	col=rep("red",length(keywords)), col.hist="grey")
{
	n<-length(keywords)
	files<-grep("\\.psmc$",dir(),value=TRUE)
	dev.new(height=height,width=width)
	for(i in 1:n)
	{
		subfiles<-grep(paste("^",keywords[i],sep=""),files,value=TRUE)
		n.sub<-length(subfiles)
		plot(1,1,
			ylim=ylim,xlim=xlim,
			log="x",type="n",
			ylab="Ne",xlab="Generations ago")
		for(i.sub in 1:n.sub)
		{
			lines(psmc.result(subfiles[i.sub],i.iteration,mu,s,g),
				type="l",col=col[i],lwd=1)
		}
		if(!is.na(command))
		{
			lines(history(command,N0,g),
				type="l",col=col.hist,lwd=2)
		}
		savePlot(filename=paste(labels[i],save.as,sep="."),type=save.as)
	}
	dev.off()
}
plotPsmc.allPops<-function(keywords, label, legend.names,
	save.as="png", height=7, width=12,
	i.iteration=25, mu=3.7e-8, s=100, g=2,
	ylim=c(0,60000), xlim=c(200,500000),
	col=rainbow(length(keywords)))
{
	n<-length(keywords)
	files<-grep("\\.psmc$",dir(),value=TRUE)
	dev.new(height=height,width=width)
	plot(1,1,ylim=ylim,xlim=xlim,log="x",type="n",main=label[1],ylab="Ne",xlab="Years ago")
	for(i in 1:n)
	{
		subfiles<-grep(paste("^",keywords[i],sep=""),files,value=TRUE)
		n.sub<-length(subfiles)
		for(i.sub in 1:n.sub)
		{
			lines(psmc.result(subfiles[i.sub],i.iteration,mu,s,g),
				type="l",col=col[i],lwd=2)
		}
	}
	legend("topright",legend=legend.names,col=col,lty=1,lwd=2)
	savePlot(filename=paste(label[1],save.as,sep="."),type=save.as)
	dev.off()
}

times<-Sys.time()
ldfile<-read.table(opt$infile,head=TRUE)
files=ldfile$file
popid=ldfile$popid
col<-rainbow(length(popid))
pop.id<-popid
psmcdata<-data.frame()
for (i in 1:length(popid)){
	data<-psmc.result(file=as.character(files[i]))
	data$Ne<-log10(data$Ne)
	data$YearsAgo<-log10(data$YearsAgo)
	data[sapply(data,is.infinite)]<-NA;
	data<-na.omit(data)
	print(c(ceiling(min(data$YearsAgo,na.rm=T)),round(max(data$YearsAgo,na.rm=T))))
	pdf(paste(popid[i],"psmc.pdf",sep="."))
	plot(1,1,ylim=c(0,round(max(data$Ne)/0.6)), xlim=c(ceiling(min(data$YearsAgo)),round(max(data$YearsAgo))),type="n",ylab="Ne(log 10)",xlab=paste("Generations ago(g=",opt$year,")",sep=""),xaxt="n")
	lines(x=data$YearsAgo,y=data$Ne,type="l",col=col[i],lwd=2)
	axis(1,at=seq(ceiling(min(data$YearsAgo)),round(max(data$YearsAgo)),1),label=paste("10","^",seq(ceiling(min(data$YearsAgo)),round(max(data$YearsAgo)),1)))
	dev.off()
	png(paste(popid[i],"psmc.png",sep="."))
	plot(1,1,ylim=c(0,round(max(data$Ne)/0.6)), xlim=c(ceiling(min(data$YearsAgo)),round(max(data$YearsAgo))),type="n",ylab="Ne(log 10)",xlab=paste("Generations ago(g=",opt$year,")",sep=""),xaxt="n")
	lines(x=data$YearsAgo,y=data$Ne,type="l",col=col[i],lwd=2)
	axis(1,at=seq(ceiling(min(data$YearsAgo)),round(max(data$YearsAgo)),1),label=paste("10","^",seq(ceiling(min(data$YearsAgo)),round(max(data$YearsAgo)),1)))
	dev.off()
	psmcdata<-rbind(psmcdata,data.frame(year=data$YearsAgo,ne=data$Ne,popid=popid[i]))
	data$YearsAgo=10^data$YearsAgo
	colnames(data)<-c(paste("Generation","x",opt$year,sep=" "),"Ne(log 10)")
	write.table(file=paste(popid[i],"psmc.result",sep="."),data,row.names=F,sep="\t")
}
pdf(paste("pop-all","psmc.pdf",sep="."))
print("haha");
plot(1,1,ylim=c(0,round(max(psmcdata$ne)/0.6)), xlim=c(ceiling(min(psmcdata$year)),round(max(psmcdata$year))),type="n",ylab="Ne(log 10)",xlab=paste("Generations ago(g=",opt$year,")",sep=""),xaxt="n")
for (i in 1:length(popid)){
	lines(y=psmcdata$ne[psmcdata$popid==popid[i]],x=psmcdata$year[psmcdata$popid==popid[i]],type="l",col=col[i],lwd=2)
}
axis(1,at=seq(ceiling(min(psmcdata$year)),round(max(psmcdata$year)),1),label=paste("10","^",seq(ceiling(min(psmcdata$year)),round(max(psmcdata$year)),1)))
legend("topright",col=col,legend=pop.id,pch=1,cex=1)
dev.off()
png(paste("pop-all","psmc.png",sep="."))
plot(1,1,ylim=c(0,round(max(psmcdata$ne)/0.6)), xlim=c(ceiling(min(psmcdata$year)),round(max(psmcdata$year))),type="n",ylab="Ne(log 10)",xlab=paste("Generations ago(g=",opt$year,")",sep=""),xaxt="n")
for (i in 1:length(popid)){
	lines(y=psmcdata$ne[psmcdata$popid==popid[i]],x=psmcdata$year[psmcdata$popid==popid[i]],type="l",col=col[i],lwd=2)
}
axis(1,at=seq(ceiling(min(psmcdata$year)),round(max(psmcdata$year)),1),label=paste("10","^",seq(ceiling(min(psmcdata$year)),round(max(psmcdata$year)),1)))
legend("topright",col=col,legend=pop.id,pch=1,cex=1)
dev.off()

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
