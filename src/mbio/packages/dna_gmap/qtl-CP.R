library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'map','m',1,'character',
	'loc','l',1,'character',
	'trt','t',1,'character',
	'out','o',1,'character',
	'num','n',1,'character',
	'method','e',1,'character',
	'pvalue','p',1,'character',
	'lod','d',1,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4)
opt = getopt(spec)
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript Rqtl_CP.r --map --loc --trt --out --num
	
Usage:
	--map	map file
	--loc	loc file
	--trt	trt file
	--out	out dir
	--num	pm number
	--method	Locating method
	--pvalue	select pvalue
	--lod	threshold lod value
	--help		usage
\n")
	q(status=1);
}
times<-Sys.time()
library('qtl');
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$map) ) { print_usage(spec) }
if ( is.null(opt$loc) ) { print_usage(spec) }
if ( is.null(opt$trt) ) { print_usage(spec) }
if ( is.null(opt$num) ) { opt$num=1000; }
if ( is.null(opt$out) ) { opt$out="./";}
if(!dir.exists(opt$out)){dir.create(opt$out)}
if (is.null(opt$method)){opt$method="cim"}

d<-read.cross(mapfile=opt$map,genfile=opt$loc,phefile=opt$trt,format="mapqtl",crosstype="4way")
setwd(opt$out);
d<-jittermap(d)
d<-sim.geno(d)
d<-calc.genoprob(d)
phe.name<-colnames(d$pheno)
nphe=length(phe.name[phe.name!= "Genotype" & phe.name !="sampleID"]);
nrow=4;
ncol=ceiling(nphe/4)
if(ncol==1){nrow=nphe}
ncol=ceiling(sqrt(nphe));
nrow=ncol;
pdf("pheno.pdf",width=30*ncol,height=40*nrow)
par(mfrow=c(ncol,nrow))
for (i in 1:length(phe.name)){
	if(phe.name[i] == "Genotype" | phe.name[i]=="sampleID"){next;}
	plotPheno(d,pheno.col=phe.name[i])
}
dev.off()
png("pheno.png")
par(mfrow=c(ncol,nrow))
for (i in 1:length(phe.name)){
	if(phe.name[i] == "Genotype" | phe.name[i]=="sampleID"){next;}
	plotPheno(d,pheno.col=phe.name[i])
}
dev.off()
chr=chrnames(d)
opt$num=as.numeric(opt$num)
opt$pvalue=as.numeric(opt$pvalue)

for(i in 1:length(phe.name)){
	if(phe.name[i] == "Genotype" | phe.name[i]=="sampleID"){next;}
	print(paste("trait",phe.name[i],sep="\t"))
	if(opt$method == "cim"){
		scan<-cim(d,pheno.col=phe.name[i]);
		scan.pm<-cim(d,pheno.col=phe.name[i],n.perm=opt$num);
	}else{
		scan<-scanone(d,pheno.col=phe.name[i]);
		scan.pm<-scanone(d,pheno.col=phe.name[i],n.perm=opt$num);
	}
	markerid<-find.marker(d,chr=scan$chr,pos=scan$pos)
	outd<-data.frame(markerid=markerid,chr=scan$chr,pos=scan$pos,lod=scan$lod);
	write.table(file=paste(phe.name[i],".scan.csv",sep=""),sep="\t",outd,row.names=FALSE)
	write.table(file=paste(phe.name[i],".pm.csv",sep=""),sep="\t",scan.pm);
	if(!is.null(opt$pvalue)){
		pm.result<-summary(scan.pm,alpha=opt$pvalue)
		write.table(file=paste(phe.name[i],".pm.summary.csv",sep=""),sep="\t",pm.result)
		scan.result<-summary(scan,format="tabByCol",threshold=pm.result,drop=1)
		threshold=pm.result
	}else{
		if(!is.null(opt$lod)){
			write.table(file=paste(phe.name[i],".pm.summary.csv",sep=""),sep="\t",pm.result)
			scan.result<-summary(scan,format="tabByCol",threshold=(opt$lod),drop=1)
			threshold=opt$lod
		}else{
			pm.result<-summary(scan.pm,alpha=0.05)
			write.table(file=paste(phe.name[i],".pm.summary.csv",sep=""),sep="\t",pm.result)
			scan.result<-summary(scan,format="tabByCol",threshold=pm.result[1,1],drop=1)
			threshold=pm.result
		}
	}
	pdf(file=paste(phe.name[i],".scan.pdf",sep=""))
	plot(scan)
	abline(h=pm.result,col=rainbow(length(pm.result)))
	dev.off()
	png(file=paste(phe.name[i],".scan.png",sep=""))
	plot(scan)
	abline(h=pm.result,col=rainbow(length(pm.result)))
	dev.off()
	if(length(rownames(scan.result$lod)) < 1){
		next;
	}
	qdata<-NULL
	n=0;
	for (j in chr){
		subd=which(outd$chr==j & outd$lod > threshold[1])
		print(paste(j,length(subd),sep="\t"))
		if(length(subd)==0){next;}
		start=subd[1]
		end=subd[1]
		if(length(subd) == 1){
			if (!is.null(qdata)){
				qdata<-rbind(qdata,data.frame(chr=j,n=n,pos=outd$pos[start:end][which.max(outd$lod[start:end])],lod=max(outd$lod[start:end]),start=outd$pos[start],end=outd$pos[end]))
			}else{
				qdata<-data.frame(chr=j,n=n,pos=outd$pos[start:end][which.max(outd$lod[start:end])],lod=max(outd$lod[start:end]),start=outd$pos[start],end=outd$pos[end])
			}
			next;
		}
		for(k in c(2:length(subd))){
			if(subd[k]-subd[k-1] < 2){
				if(subd[k-1] < start){start=subd[k-1]}
				if(subd[k] > end){end=subd[k]}
			}else{	
				print(paste(start,end))
				if (!is.null(qdata)){
					qdata<-rbind(qdata,data.frame(chr=j,n=n,pos=outd$pos[start:end][which.max(outd$lod[start:end])],lod=max(outd$lod[start:end]),start=outd$pos[start],end=outd$pos[end]))
				}else{
					qdata<-data.frame(chr=j,n=n,pos=outd$pos[start:end][which.max(outd$lod[start:end])],lod=max(outd$lod[start:end]),start=outd$pos[start],end=outd$pos[end])
				}
				n=n+1
				start=subd[k]
				end=subd[k]
			}
		}
		if(end != subd[1]){
			n=n+1;
			if (!is.null(qdata)){
				qdata<-rbind(qdata,data.frame(chr=j,n=n,pos=outd$pos[start:end][which.max(outd$lod[start:end])],lod=max(outd$lod[start:end]),start=outd$pos[start],end=outd$pos[end]))
			}else{
				qdata<-data.frame(chr=j,n=n,pos=outd$pos[start:end][which.max(outd$lod[start:end])],lod=max(outd$lod[start:end]),start=outd$pos[start],end=outd$pos[end])
			}
		}
	}
	qtlname=paste(phe.name[i],c(1:length(qdata$n)),sep="-")
	qtl<-makeqtl(d,chr=qdata$chr,pos=qdata$pos,qtl.name=qtlname)
	fitqtl<-fitqtl(cross=d,qtl=qtl,get.est=TRUE,pheno.col=i)
	markerid<-find.marker(d,chr=qtl$chr,pos=qtl$pos)
	var<-fitqtl$result.drop[,"%var"]
	if (length(qtl$name) == 1){var<-fitqtl$result.full["Model","%var"]}
	data<-data.frame(marker=markerid,chr=qdata$chr,pos=qdata$pos,lod=qdata$lod,var=var,pm1=pm.result[1],pm2=pm.result[2],start=qdata$start,end=qdata$end)
	for(j in 1:length(qtlname)){
		data$mark1[j]=find.marker(d,chr=qtl$chr[j],data$start[j])
		data$mark2[j]=find.marker(d,chr=qtl$chr[j],data$end[j])
		pdf(paste(phe.name[i],".",qtl$name[j],".PXG.pdf",sep=""))
		plotPXG(d,data$marker[j],pheno.col=i)
		dev.off()
		png(paste(phe.name[i],".",qtl$name[i],".PXG.png",sep=""))
		plotPXG(d,data$marker[j],pheno.col=i)
		dev.off()
	}
	write.table(file=paste(phe.name[i],".qtl.csv",sep=""),sep="\t",data,row.names=FALSE)
	pdf(paste(phe.name[i],".qtl.pdf",sep=""))
	plot(qtl)
	dev.off()
	png(paste(phe.name[i],".qtl.png",sep=""))
	plot(qtl)
	dev.off()
}
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
