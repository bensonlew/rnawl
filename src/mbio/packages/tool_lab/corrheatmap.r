# # # # ~/app/bioinfo/rna/miniconda2/bin/R

# options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# install.packages("psych")
# install.packages("RColorBrewer")
# install.packages("corrplot")

library('getopt');

spec = matrix(c(
	'infile','i',1,'character',
	'outdir','o',1,'character',
	'type','t',2,'character',               #"full", "lower", "upper"
	'method','m',2,'character',             #"circle", "square", "ellipse", "number", "shade", "color", "pie"
	'color','c',2,'character',				 # Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd, BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral
	'corr_method','d',2,'character',        # "pearson", "kendall", "spearman"
	'corr_adjust','a',2,'character',        #"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
	'show_significance','s',2,'character',  # 'yes','no'
	'in_significance','n',2,'character',    # 'pch' (default), 'p-value', 'blank', 'label_sig'
	'show_coefficient','e',2,'character',    # 'yes','no'
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--infile(-i)	 the input file
	--outdir(-o)	 the output directory
	--type(-t)	     Character, 'full' (default), 'upper' or 'lower', display full matrix, lower triangular or upper triangular matrix.
	--method(-m)	 Character, the visualization method of correlation matrix to be used: circle, square, ellipse, number, shade, color, pie
	--color(-c)	     A palette name from the lists below: Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd, BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral
	--corr_method(-d)	relation method: pearson, kendall, spearman
	--corr_adjust(-a)	The adjustment for multiple tests should be used : holm, hochberg, hommel, bonferroni, BH, BY, fdr, none
	--show_significance(-s)	whether show significance:  yes,no
	--in_significance(-n)	the way to show significance:  pch (default), p-value, blank, label_sig
	--show_coefficient(-e)	whether show coefficient on plot:  yes,no
	--help(-h)		usage
\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$infile)) { print_usage(spec)}
if ( is.null(opt$outdir)) { print_usage(spec)}
if ( is.null(opt$type)){ opt$type="full"}
if ( is.null(opt$method)){ opt$method="circle"}
if ( is.null(opt$color)){ opt$color="Spectral" }
if ( is.null(opt$corr_method)){ opt$corr_method="pearson"}
if ( is.null(opt$corr_adjust)){ opt$corr_adjust="holm"}
if ( is.null(opt$show_significance)){ opt$show_significance="yes"}
if ( is.null(opt$in_significance)){ opt$in_significance="pch"}
if ( is.null(opt$show_coefficient)){ opt$show_coefficient="yes"}

print(paste("infile ",opt$infile))
print(paste("outdir： ",opt$outdir))
print(paste("type ",opt$type))
print(paste("method ",opt$method))
print(paste("color ",opt$color))
print(paste("corr_method ",opt$corr_method))
print(paste("corr_adjust ",opt$corr_adjust))
print(paste("show_significance ",opt$show_significance))
print(paste("in_significance ",opt$in_significance))
print(paste("show_coefficient ",opt$show_coefficient))

data = as.data.frame((read.table(opt$infile,header=T,sep="\t",row.names=1)))

library("psych")
library("RColorBrewer")
library("corrplot")

runR = paste( "corr_matrix <- corr.test(data, method = '",opt$corr_method,"', adjust = '", opt$corr_adjust ,"')" ,sep="")
print('excute R command1:')
print(runR)
print('Finish excute R command1')
eval(parse(text = runR ) )

pdf(file.path(opt$outdir, 'corr_heatmap.pdf'), width = 13, height = 13)
runR = paste( "p <- corrplot(corr_matrix$r , method = '", opt$method ,"', type = '", opt$type ,"', col = brewer.pal(n=8, name = '", opt$color ,"'), ",sep="")
if (opt$show_significance == "yes"){
	runR = paste(runR,"p.mat = corr_matrix$p, insig = '", opt$in_significance ,"',pch.cex = 0.9,pch.col = 'grey20',sig.level = c(0.001, 0.01, 0.05), ",sep="")
}
if (opt$show_coefficient == "yes"){
	runR = paste(runR,"addCoef.col = 'black', number.font= 5", sep="")
}
runR = paste(runR,")", sep="")
print('excute R command2:')
print(runR)
print('Finish excute R command2')
eval(parse(text = runR ))

dev.off()
# library("ggplot2")
# ggsave(filename = file.path(opt$outdir, 'corr_heatmap.pdf'), plot = p, width = 13, height = 13)
