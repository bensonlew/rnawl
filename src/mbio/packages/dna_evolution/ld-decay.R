#!/usr/bin/env Rscript
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'infile','i',0,'character',
	'outfile','o',0,'character',
	'help','m',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	q(status=1);
}
if ( !is.null(opt$help)) { print_usage(spec) }
if ( is.null(opt$infile)) { print_usage(spec)}
if ( is.null(opt$outfile)){ print_usage(spec) }

times<-Sys.time()
ldfile<-read.table(opt$infile,head=TRUE)
files=ldfile$file
popid=ldfile$popid
num=ldfile$sampleNum
col<-rainbow(length(popid))
#pop.id<-paste("pop",c(1:length(popid)))
newLD<-NULL
drawdis<-c(0:300000)
out<-NULL
for (i in 1:length(popid)){
	ld=read.table(file=as.character(files[i]),head=TRUE,comment.char=":");
	distance=ld$X.Dist
	R2=ld$Mean_r.2
	n=num[i]
	HW.st<-c(C=0.0001)
	HW.nonlinear<-nls(R2~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=1000))
	tt<-summary(HW.nonlinear)
	new.rho<-tt$parameters[1]
	maxld<-max(ld$Mean_r.2)
	decay05<-distance[which.min(abs(ld$Mean_r.2-maxld/2))]
	decay01<-distance[which.min(abs(ld$Mean_r.2-0.1))]
	distance=drawdis
	fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
	fpoints[1]=max(R2)
	if (is.null(newLD)){
		newLD=data.frame(distance/1000,fpoints)
	}else{
		newLD=cbind(newLD,fpoints)
	}
	maxld.1<-max(fpoints)
	decay05.1<-fpoints[which.min(abs(fpoints-fpoints-1/2))]
	decay01.1<-fpoints[which.min(abs(fpoints-0.1))]

	out<-rbind(out,data.frame(popid[i],decay05,decay01,maxld,decay05.1,decay01.1,maxld.1))
}
colnames(newLD)=c("distance",popid)
png(paste(opt$outfile,"png",sep="."));
for (i in 2:(length(popid)+1)){
	if (i!=2){
		lines(x=newLD$distance,y=newLD[,i],col=col[i-1],lwd=1)
	}else{
		plot(x=newLD$distance,y=newLD[,i],type="l",lwd=1,ylim=c(0,1),col=col[i-1],main="LD decay",ylab="R^2",xlab="distance(kb)")
		#lines(x=newLD$distance,y=rev(newLD[,i]),col=col[i-1],lwd=1)
	}
}
legend("topright",col=col,legend=popid,pch=1,cex=1)
dev.off();
pdf(paste(opt$outfile,"pdf",sep="."));
for (i in 2:(length(popid)+1)){
	if (i!=2){
		lines(x=newLD$distance,y=newLD[,i],col=col[i-1],lwd=1)
	}else{
		plot(x=newLD$distance,y=newLD[,i],type="l",lwd=1,ylim=c(0,1),col=col[i-1],main="LD decay",ylab="R^2",xlab="distance(kb)")
		#lines(x=newLD$distance,y=newLD[,i],col=col[i],lwd=1)
	}
}
legend("topright",col=col,legend=popid,pch=1,cex=1)
dev.off();

write.table(file=paste(opt$outfile,"table",sep="."),out)
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
