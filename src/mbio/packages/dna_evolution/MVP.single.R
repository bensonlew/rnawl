library('getopt');

spec = matrix(c(
	'hapmap','h',0,'character',
	'trait','t',0,'character',
	'output','o',0,'character',
	'help','c',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript indel_len.r --i  --o  
	
Usage:
	--hapmap	the input hapmap file
	--trait	the trait file 
	--output	the output dir
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$trait)) { print_usage(spec)}
if ( is.null(opt$hapmap)){ print_usage(spec) }
if ( is.null(opt$output)){ print_usage(spec) }

library(MVP)
setwd(opt$output);
MVP.Data(fileHMP=opt$hapmap,
         filePhe=opt$trait,
         sep.hmp="\t",
         sep.phe="\t",
         SNP.effect="Add",
         fileKin=FALSE,
         filePC=FALSE,
         out="mvp.hmp",
         priority="memory",
)
genotype <- attach.big.matrix("mvp.hmp.geno.desc")
phenotype <- read.table("mvp.hmp.phe",head=TRUE)
map <- read.table("mvp.hmp.map" , head = TRUE)
#MVP.Hist(phe=phenotype[,c(1,i)], file="png", breakNum=30, dpi=300)
imMVP <- MVP(
		phe=phenotype,
		geno=genotype,
		map=map,
		nPC.GLM=5,
		nPC.MLM=3,
		nPC.FarmCPU=3,
		perc=1,
		priority="speed",
		ncpus=10,
		vc.method="EMMA",
		maxLoop=10,
		method.bin="FaST-LMM",
		threshold=0.05,
		file="pdf",
		method=c("GLM", "MLM", "FarmCPU")
)
#MVP.Report(imMVP, plot.type="m", LOG10=TRUE, ylim=NULL, threshold=c(1e-6,1e-4),threshold.lty=c(1,2),col=c("grey60","grey30"), threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,chr.den.col=c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3"),bin.size=1e6,signal.col=c("red","green"),signal.cex=c(1,1),signal.pch=c(19,19),file="jpg",memo="",dpi=300)

