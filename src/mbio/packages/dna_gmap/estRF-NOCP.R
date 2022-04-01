library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'mark','m',1,'character',
	'trt','t',1,'character',
	'out','o',1,'character',
	'num','n',1,'character',
	'method','d',1,'character',
	'pvalue','q',1,'character',
	'lod','l',1,'character',
	'pop','p',1,'character',
	'bc','b',1,'character',
	'f','f',1,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4)
opt = getopt(spec)
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript Rqtl-NOCP.r --mark  --out --num --pop --method
	
Usage:
	--mark	map file
	--trt	trt file
	--pop	pop type
	--out	out dir
	--num	pm number
	--method	method
	--pvalue	pvalue
	--lod	lod
	--bc	bc gen for bcsft
	--f		f gen for bcsft
	--help		usage
\n")
	q(status=1);
}
times<-Sys.time()
library('qtl');
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$mark) ) { print_usage(spec) }
if ( is.null(opt$trt) ) { print_usage(spec) }
if ( is.null(opt$pop) ) { print_usage(spec) }
if ( is.null(opt$num) ) { opt$num=1000; }
if ( is.null(opt$out) ) { opt$out="./";}
if(opt$pop =="bcsft" & is.null(opt$bc) & is.null(opt$f)){print_usage(spec)}
if(opt$pop =="bcsft") {	
	d<-read.cross(file=opt$mark,phefile=opt$trt,format="csvsr",crosstype=opt$pop,na.string="NaN",BC.gen=opt$bc,F.gen=opt$f)
}else{
	d<-read.cross(file=opt$mark,phefile=opt$trt,format="csvsr",crosstype=opt$pop,na.string="NaN")
}
setwd(opt$out);
d<-est.rf(d)
write.table(file="pop.rf",d$rf)
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
