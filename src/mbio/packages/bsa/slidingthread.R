#!/share/nas2/genome/biosoft/R/2.15.1/lib64/R/bin/Rscript 
# load library
times<-Sys.time()
library('getopt');
options(bitmapType='cairo')
#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
	'help' , 'h', 0, "logical",
	'infile' , 'i', 1, "character",
	'outfile' , 'o', 1, "character",
	'thresold','t',1,"character",
	'col' , 'c', 1, "integer"
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
# define usage function
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
Options: 
	--help		NULL 		get this help
	--infile 	character 	the input file [forced]
	--outfile 	character 	the filename for output graph [forced]
	--col		character	the col to thread default 4
	--thresold	character	the threhold q value eg \"0.95,0.99,0.999\"
	\n")
	q(status=1);
}
# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }
if (  is.null(opt$col) ) { opt$col=4 }
# check non-null args
if ( is.null(opt$infile) )	{ print_usage(spec) }
if ( is.null(opt$outfile) )	{ print_usage(spec) }
if ( is.null(opt$cutoff))	{opt$cutoff="0.95,0.99,0.999"}
if ( opt$col == 4){
	idxsmooth=read.table(opt$infile,header=TRUE,col.names=c("chr","start","end","idx","snpnum"));
}else{
	idxsmooth=read.table(opt$infile,header=TRUE,col.names=c("chr","start","end","pdx","mdx","idx","snpnum"));
}
cutoff=apply(as.matrix(unlist(strsplit(opt$cutoff,split=","))),2,as.numeric);
# output quantile index
if(file.exists(paste(opt$outfile,".quantile",sep=""))){
	file.remove(paste(opt$outfile,".quantile",sep=""))
}
if(file.exists(paste(opt$outfile,".result",sep=""))){
	file.remove(paste(opt$outfile,".result",sep=""))
}
for (i in 1:9999){
	quant=quantile(idxsmooth$idx,i/10000);
	nlen=length(idxsmooth$idx)*(1-i/10000);
	newdata <- cbind(i/10000, quant,nlen);
	write.table(newdata,file=paste(opt$outfile,".quantile",sep=""),append = TRUE, row=F,col=F,sep ="\t",quote=F);
}
data<-cbind(cutoff,quantile(idxsmooth$idx,cutoff))
write.table(data,file=paste(opt$outfile,".result",sep=""),append = TRUE, row=F,col=F,sep ="\t",quote=F);	
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)


