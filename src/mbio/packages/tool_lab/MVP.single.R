library('getopt');

spec = matrix(c(
	'plink','p',0,'character',
	'output','o',0,'character',
	'help','c',0,'logical',
	# 'low','l',0,'character',
	# 'middle','m',0,'character',
	# 'high','h',0,'character',
	'wsize','w',0,'character'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript indel_len.r --i  --o  
	
Usage:
	--plink	the input plink file
	--output	the output dir
	--low       low colour
	--middle    middle colour
	--high      high colour
	--wsize	    window size
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$plink)){ print_usage(spec) }
if ( is.null(opt$output)){ print_usage(spec) }
# if ( is.null(opt$low)){opt$low = "yellow"}
# if ( is.null(opt$high)){opt$high="red"}
# if (is.null(opt$middle)){opt$middle="green"}
if(is.null(opt$wsize)){opt$wsize=1e6}
library(rMVP)
setwd(opt$output);

MVP.Data(fileBed=opt$plink,
         fileKin=FALSE,
         filePC=FALSE,
         out="mvp.plink",         
         priority="memory",
         )

map <- read.table("mvp.plink.geno.map" , head = TRUE)
MVP.Report(map[,1:3], plot.type="d",  bin.size=as.numeric(opt$wsize), file.type="pdf")
MVP.Report(map[,1:3], plot.type="d",  bin.size=as.numeric(opt$wsize), file.type="jpg")
