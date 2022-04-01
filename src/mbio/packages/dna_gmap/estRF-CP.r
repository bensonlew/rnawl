library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'map','m',1,'character',
	'loc','l',1,'character',
	'trt','t',1,'character',
	'out','o',1,'character',
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
if ( is.null(opt$out) ) { opt$out="./";}
if(!dir.exists(opt$out)){dir.create(opt$out)}

d<-read.cross(mapfile=opt$map,genfile=opt$loc,phefile=opt$trt,format="mapqtl",crosstype="4way")
setwd(opt$out);
d<-est.rf(d)
#cat(paste(opt$out, "pop.rf", sep = "/"))
write.table(file= "pop.rf", d$rf)

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
