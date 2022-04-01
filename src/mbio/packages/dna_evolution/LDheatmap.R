#!/usr/bin/env Rscript
library('getopt');
library(ggplot2)
library(grid)

options(bitmapType='cairo')
spec = matrix(c(
	'input','i',0,'character',
	'output','o',0,'character',
	'help','m',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--input	the input  file
	--output	the out file
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$input)) { print_usage(spec)}
if ( is.null(opt$output)){ print_usage(spec) }
times<-Sys.time()
#setwd(opt$output)
#############################################
library(genetics)
library(LDheatmap)
data<-read.table(opt$input,header=T)
genematrix<-data[-1,]
genedis<-as.numeric(t(data[1,]))
num<-ncol(genematrix)

for (i in 1:num){
	genematrix[,i]<-as.genotype(genematrix[,i])
}

pdf(paste(opt$output,"LDheatmap.pdf",sep="."))
MyHeatmap<-LDheatmap(genematrix,genetic.distances=as.vector(t(genedis)),color = c("#08306B","#084D96","#1B69AF","#3787C0","#58A1CE","#81BADA","#ABCFE5","#CBDEF0","#E0ECF7","#F7FBFF"),flip=TRUE,title="Pairwise LD",add.map=TRUE, add.key=TRUE,distances="physical")
dev.off()

#flippedHeatmap<-LDheatmap(MyHeatmap,flip=TRUE) ##翻转图片使其成倒三角
#LDheatmap(MyHeatmap, SNP.name = c("**", "**")) ## 显示**和**标记名称
#grid.edit(gPath("ldheatmap", "heatMap", "title"), gp = gpar(col = "red")) ##修改标题颜色，可加大小等
#grid.edit(gPath("ldheatmap", "geneMap","SNPnames"), gp = #gpar(cex=1.5)) ##修改SNP标签大小
#grid.edit(gPath("ldheatmap", "heatMap", "heatmap"), gp = gpar(col = "white",lwd = 2))  ##使用白色线条间隔LD色块


png(paste(opt$output,"LDheatmap.png",sep="."))
MyHeatmap<-LDheatmap(genematrix,genetic.distances=as.vector(t(genedis)),color = c("#08306B","#084D96","#1B69AF","#3787C0","#58A1CE","#81BADA","#ABCFE5","#CBDEF0","#E0ECF7","#F7FBFF"),flip=TRUE,title="Pairwise LD",add.map=TRUE, add.key=TRUE,distances="physical")
dev.off()

#############################################
escaptime=Sys.time()-times;

print("Done!")
print(escaptime)
