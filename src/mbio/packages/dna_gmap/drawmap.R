library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'mark','m',1,'character',
	'out','o',1,'character',
	'pop','p',1,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4)
opt = getopt(spec)
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript emap.R --mark --out --pop
	
Usage:
	--mark	map file
	--out	out dir
	--pop	pop type
	--help		usage
\n")
	q(status=1);
}
times<-Sys.time()
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$mark) ) { print_usage(spec) }
if ( is.null(opt$pop) ) { print_usage(spec) }

if ( is.null(opt$out) ) { opt$out="./";}
if(!dir.exists(opt$out)){dir.create(opt$out)}
library('qtl');
library('ASMap');
library('Cairo')


opt$pop=tolower(opt$pop)
if(opt$pop == "cp"){
	d<-read.cross(genfile=paste(opt$mark,"sexAver.loc",sep="."),phefile=paste(opt$mark,"trt",sep="."),mapfile=paste(opt$mark,"sexAver.map",sep="."),format="mapqtl")
	d1<-read.cross(genfile=paste(opt$mark,"sexAver.loc",sep="."),phefile=paste(opt$mark,"trt",sep="."),mapfile=paste(opt$mark,"male.map",sep="."),format="mapqtl")
	d2<-read.cross(genfile=paste(opt$mark,"sexAver.loc",sep="."),phefile=paste(opt$mark,"trt",sep="."),mapfile=paste(opt$mark,"female.map",sep="."),format="mapqtl")
	setwd(opt$out)
	d<-jittermap(d)
	d1<-jittermap(d1)
	d2<-jittermap(d2)
	# pdf("total.sexAver.pdf");
	# plotMap(d,shift=TRUE,alternate.chrid=TRUE)
	# dev.off()
	# png("total.sexAver.png");
	# plotMap(d,shift=TRUE,alternate.chrid=TRUE)
	# dev.off()
	# pdf("total.male.pdf");
	# plotMap(d1,shift=TRUE,alternate.chrid=TRUE)
	# dev.off()
	# png("total.male.png");
	# plotMap(d1,shift=TRUE,alternate.chrid=TRUE)
	# dev.off()
	# pdf("total.female.pdf");
	# plotMap(d2,shift=TRUE,alternate.chrid=TRUE)
	# dev.off()
	# png("total.female.png");
	# plotMap(d2,shift=TRUE,alternate.chrid=TRUE)
	# dev.off()
	chrname<-chrnames(d);
	for(i in c(1:length(chrname))){
		pdf(paste(chrname[i],".heatMap.sexaver.pdf",sep=""))
		plotRF(d,chr=i,lmax=50,"both")
		dev.off()
		CairoPNG(paste(chrname[i],".heatMap.sexaver.png",sep=""))
		plotRF(d,chr=i,lmax=50,"both")
		dev.off()
	}
	chrname<-chrnames(d1);
	for(i in c(1:length(chrname))){
		pdf(paste(chrname[i],".heatMap.male.pdf",sep=""))
		plotRF(d1,chr=i,lmax=50,"both")
		dev.off()
		CairoPNG(paste(chrname[i],".heatMap.male.png",sep=""))
		plotRF(d1,chr=i,lmax=50,"both")
		dev.off()
	}
	chrname<-chrnames(d2);
	for(i in c(1:length(chrname))){
		pdf(paste(chrname[i],".heatMap.female.pdf",sep=""))
		plotRF(d2,chr=i,lmax=50,"both")
		dev.off()
		CairoPNG(paste(chrname[i],".heatMap.female.png",sep=""))
		plotRF(d2,chr=i,lmax=50,"both")
		dev.off()
	}
}else{
	d<-read.cross(file=paste(opt$mark,"csv",sep="."),format="csvr",crosstype=opt$pop)
	setwd(opt$out)
	d<-jittermap(d)
	d<-est.rf(d)
	# pdf("total.lg.pdf");
	# plotMap(d,shift=TRUE,alternate.chrid=TRUE)
	# dev.off()
	# png("total.lg.png");
	# plotMap(d,shift=TRUE,alternate.chrid=TRUE)
	# dev.off()
	chrname<-chrnames(d);
	for(i in c(1:length(chrname))){
		pdf(paste(i,".heatMap.pdf",sep=""))
		plotRF(d,chr=i,what='both')
		dev.off()
		CairoPNG(paste(i,".heatMap.png",sep=""))
		plotRF(d,chr=i,what='both')
		dev.off()
	}
}

escaptime=Sys.time()-times;
print("Done!\n")
print(escaptime)
