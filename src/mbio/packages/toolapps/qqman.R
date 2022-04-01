#!/usr/bin/env Rscript
#vignette('qqman')
library('qqman');
library('getopt');
spec = matrix(c(
				'snp','s',1,'character',
				'suggestiveline','l',2,'numeric',
				'secondline','y',2,'character',
				'genomewideline','g',2,'numeric',
				'highlight','t',2,'character',
				'col','c',2,'character',
				'help', 'h', 0, 'logical',
				'picone','p',2,'numeric',
				'chrlabs','b',2,'character',
				'outdir','o',2,'character'
				), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
cat("
		Usage example: 
		Rscript DESeq2.R -m counts.txt -g group.txt -c condition.txt -o outdir 
				Usage:
                   -s   --snp                    snp_info_file with clumus CHR,BP,P
                   -l   --suggestiveline         defalse 1e-5 because y=-log10(1e-5)
				   -y   --secondline             draw the secondline or not values(yes,no)
				   -g	--genomewideline         defalse 5e-8 because y=-log10(5e-8)
				   -t	--highlight	             highlight snp_id need four columns
				   -c   --col                    two colour you want
				   -b   --chrlabs                axis_x_labs if your chr has been change
				   -p   --picone                 Show one chromosome
				   -o   --outdir                 output_dir
				   -h   --help	                 usage\n")
	q(status=1);}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$snp) ) { print_usage(spec) }
if ( is.null(opt$suggestiveline) ) { opt$suggestiveline = -log10(1e-5) }else if(opt$suggestiveline == 0){opt$suggestiveline=FALSE}else{opt$suggestiveline = -log10(opt$suggestiveline)}
if ( is.null(opt$secondline) ) { opt$secondline = "no" }else if(opt$secondline == "yes"){opt$secondline = "yes"}
if ( is.null(opt$genomewideline) ) { opt$genomewideline = -log10(5e-8)}else{opt$genomewideline=-log10(opt$genomewideline)}
if ( is.null(opt$outdir) ) { opt$outdir="./" }
if ( is.null(opt$col)){col_list=c("#A4D3EE","#DDA0DD")}else{col_list=c(strsplit(opt$col,",")[[1]])}
if(!dir.exists(opt$outdir)){dir.create(opt$outdir)}
SNP_info = read.table(opt$snp,header=T,sep="\t")
if(opt$secondline == "yes"){opt$genomewideline=opt$genomewideline}else{opt$genomewideline=FALSE}
if(length(SNP_info)==3){
colnames(SNP_info)=c("CHR","BP","P")
}else if(length(SNP_info)==4){
colnames(SNP_info)=c("SNP","CHR","BP","P")
}
if(!is.null(opt$chrlabs)){chrlabs_list = as.character(c(strsplit(opt$chrlabs,",")[[1]]))
}

if(is.null(opt$picone)){
  if(is.null(opt$chrlabs)){
	if(is.null(opt$highlight)){
    png(paste0(opt$outdir,"/qqman.png"))	
	manhattan(SNP_info,col=col_list,suggestiveline=opt$suggestiveline,genomewideline=opt$genomewideline)

	dev.off()
	pdf(paste0(opt$outdir,"/qqman.pdf"))
	manhattan(SNP_info,col=col_list,suggestiveline=opt$suggestiveline,genomewideline=opt$genomewideline)
	dev.off()
	}else{
	snpsOfInterest = c(strsplit(opt$highlight,",")[[1]])
	png(paste0(opt$outdir,"/qqman.png"))
	manhattan(SNP_info,highlight = snpsOfInterest,col=col_list,suggestiveline=opt$suggestiveline,genomewideline=opt$genomewideline)
	dev.off()
	pdf(paste0(opt$outdir,"/qqman.pdf"))
	manhattan(SNP_info,highlight = snpsOfInterest,col=col_list,suggestiveline=opt$suggestiveline,genomewideline=opt$genomewideline)
	dev.off()
	}
    }else{
  
       if(is.null(opt$highlight)){
		   png(paste0(opt$outdir,"/qqman.png"))
		   manhattan(SNP_info,col=col_list,suggestiveline=opt$suggestiveline,genomewideline=opt$genomewideline,chrlabs=chrlabs_list)
		   dev.off()
		   pdf(paste0(opt$outdir,"/qqman.pdf"))
		   manhattan(SNP_info,col=col_list,suggestiveline=opt$suggestiveline,genomewideline=opt$genomewideline,chrlabs=chrlabs_list)
		   dev.off()
	   }else{
		   snpsOfInterest = c(strsplit(opt$highlight,",")[[1]])
		   png(paste0(opt$outdir,"/qqman.png"))
		   manhattan(SNP_info,highlight = snpsOfInterest,col=col_list,suggestiveline=opt$suggestiveline,genomewideline=opt$genomewideline,chrlabs=chrlabs_list)
		   dev.off()
		   pdf(paste0(opt$outdir,"/qqman.pdf"))
		   manhattan(SNP_info,highlight = snpsOfInterest,col=col_list,suggestiveline=opt$suggestiveline,genomewideline=opt$genomewideline,chrlabs=chrlabs_list)
		   dev.off()
	   }
  }
######################################
}else{
   if(is.null(opt$chrlabs)){
   png(paste0(opt$outdir,"/qqman.png"))
   manhattan(subset(SNP_info,CHR== opt$picone),col=col_list[1],suggestiveline=opt$suggestiveline,genomewideline=opt$genomewideline)
   dev.off()
   pdf(paste0(opt$outdir,"/qqman.pdf"))
   manhattan(subset(SNP_info,CHR==opt$picone),col=col_list[1],suggestiveline=opt$suggestiveline,genomewideline=opt$genomewideline)
   dev.off()}else{
   png(paste0(opt$outdir,"/qqman.png"))
   manhattan(subset(SNP_info,CHR== opt$picone),col=col_list[1],suggestiveline=opt$suggestiveline,genomewideline=opt$genomewideline,chrlabs=chrlabs_list)
   dev.off()
   pdf(paste0(opt$outdir,"/qqman.pdf"))
   manhattan(subset(SNP_info,CHR==opt$picone),col=col_list[1],suggestiveline=opt$suggestiveline,genomewideline=opt$genomewideline,chrlabs=chrlabs_list)
   dev.off()
   }
}



