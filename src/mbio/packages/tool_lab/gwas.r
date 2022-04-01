library('getopt');

spec = matrix(c(
	'plink','h',0,'character',
	'trait','t',0,'character',
	'output','o',0,'character',
	'help','c',0,'logical',
	'method','m',0,'character',
	'threshold','p',0,'character',
	'chrom', 'r', 0, 'character'
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
	--trait	the trait file 
	--output	the output dir
	--method	the method \"GLM\",\"FarmCPU\",\"MLM\"
	--threshold the threshold default 0.05
	--chrom     the chromosome map file
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$trait)) { print_usage(spec)}
if ( is.null(opt$plink)){ print_usage(spec) }
if ( is.null(opt$output)){ print_usage(spec) }
if ( is.null(opt$method)){opt$method = "MLM"}
if ( is.null(opt$threshold)){opt$threshold=0.05}
library(rMVP)
setwd(opt$output);

MVP.Data(fileBed=opt$plink,
         filePhe=opt$trait,
         fileKin=FALSE,
         filePC=FALSE,
         out="mvp.plink",         
         priority="memory",
         )

genotype <- attach.big.matrix("mvp.plink.geno.desc")
phenotype <- read.table("mvp.plink.phe",head=TRUE)
map <- read.table("mvp.plink.geno.map" , head = TRUE)
chrom<-read.table(opt$chrom)
#MVP.Hist(phe=phenotype[,c(1,i)], file="png", breakNum=30, dpi=300)

for (i in c(2:ncol(phenotype)-1)){

imMVP <- MVP(phe=phenotype[,c(1,i+1)],
    geno=genotype,
	map=map,
	nPC.GLM=5,
	nPC.MLM=3,
	nPC.FarmCPU=3,
	priority="speed",
	ncpus=10,
	vc.method="EMMA",
	maxLoop=10,
	method.bin="FaST-LMM",
	p.threshold=opt$threshold,
	file.output=FALSE,
	method=c(opt$method)
    )

    result <- cbind(imMVP$map[,1:3], imMVP[[paste0(tolower(opt$method), '.results')]][,3])
    sample<-colnames(result)[4]<-colnames(imMVP[[paste0(tolower(opt$method), '.results')]])[3]
    order <- order(result[,1], result[,2])
    order_chr<-unique(result[order,2])

    chrom_v<-c()
    for (i in order_chr){
        chr<-chrom[chrom[[2]] == i,1]
        chrom_v<-c(chrom_v, as.character(chr))
    }

    MVP.Report(
        result,
        plot.type=c("m","q"),
        file.output=TRUE,
        file.type='pdf',
        chr.labels = chrom_v,
        chr.den.col=c("darkgreen", "yellow", "red"),
        dpi=300,
        threshold=as.numeric(opt$threshold)/nrow(genotype),
    )
    result<-cbind(imMVP$map[,1:3], imMVP[[paste0(tolower(opt$method), '.results')]][,c(1,3)])
    result[,2]<-sapply(result[,2], function(x) chrom[chrom[[2]] == x,1])
    colnames(result)<-c('snp', 'chrom', 'bp', 'effect', 'pvalue')
    write.table(result, paste0(sample,'_gwas_loci.xls'), sep='\t', quote=F, row.names=F)
}
    #lamda<-function(qvalue=NULL){
    #	z=qnorm(qvalue/2)
    #	return(round(median(z^2, na.rm = TRUE) / qchisq(0.5, 1), 3))
    #}
    #colnames(imMVP$glm.results)[1:2]=c("G.effect","G.se")
    #colnames(imMVP$mlm.results)[1:2]=c("M.effect","M.se")
    #colnames(imMVP$farmcpu.results)[1:2]=c("F.effect","F.se")
    #iMVP.res <- cbind(map,imMVP$glm.results, imMVP$mlm.results, imMVP$farmcpu.results)
    #imMVP<-select(iMVP.res,1:3|ends_with("MLM") | ends_with("GLM")|ends_with("FarmCPU"))
    #lamdas<-apply(imMVP[,3:ncol(imMVP)],2,lamda)
    #write.table(file=paste0(sample,"lamda.xls"),lamdas,sep="\t",quote="")
    #write.table(file=paste0(sample, "imMVP.merge.result"),imMVP,sep="\t",quote="",row.names=F)
    #MVP.Report(imMVP, plot.type="m", LOG10=TRUE, ylim=NULL, threshold=c(1e-6,1e-4),threshold.lty=c(1,2),col=c("grey60","grey30"), threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,chr.den.col=c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3"),bin.size=1e6,signal.col=c("red","green"),signal.cex=c(1,1),signal.pch=c(19,19),file="jpg",memo="",dpi=300)

