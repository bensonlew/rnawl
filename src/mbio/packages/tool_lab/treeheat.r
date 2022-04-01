# # # # ~/app/bioinfo/ref_rna_v3/clusterprofile_4.1/miniconda3/bin/R

# export PATH=/mnt/lustre/users/sanger-dev/app/bioinfo/ref_rna_v3/clusterprofile_4.1/miniconda3/bin/:$PATH
# options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# install.packages("treeheatr")

library('getopt');

spec = matrix(c(
	'infile','i',1,'character',
	'outdir','o',1,'character',
	'target_lab','t',2,'character',               # column name which include No NA
	'show_apart','s',2,'character',             # "heat-tree", "heat-only", "tree-only"
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
	--target_lab(-t)	     Character, column name which include No NA.
	--show_apart(-s)	 Character, Character string indicating which components of the decision tree-heatmap should be drawn. Can be ’heat-tree’, ’heat-only’ or ’tree-only’.
	--help(-h)		usage
\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$infile)) { print_usage(spec)}
if ( is.null(opt$outdir)) { print_usage(spec)}
if ( is.null(opt$target_lab)){ print_usage(spec) }
if ( is.null(opt$show_apart)){ opt$show_apart="heat-tree"}

print(paste("infile ",opt$infile))
print(paste("outdir： ",opt$outdir))
print(paste("target_lab ",opt$target_lab))
print(paste("show_apart ",opt$show_apart))

data = as.data.frame((read.table(opt$infile,header=T,sep="\t")))
good_columns_with_no_na = names(which(colSums(is.na(data)) == 0))

if (opt$target_lab %in% good_columns_with_no_na){
	library("treeheatr")
	pdf(file.path(opt$outdir, 'treeheat.pdf'), width = 12, height = 9, onefile=FALSE)
	heat_tree(data, target_lab = opt$target_lab,show = opt$show_apart,cont_legend = TRUE,cate_legend = TRUE)
}
if (opt$target_lab %in% good_columns_with_no_na){
	dev.off()
}


# tryCatch({
# 	pdf(file.path("C:/Users/xi.xu/Desktop/热图决策树--小工具", 'test14.pdf'), width = 12, height = 6, onefile=FALSE)
# 	heat_tree(xx_df, target_lab = 'sex',show = 'heat-tree')
# 	dev.off()
# 	}, warning = function(w){
# 		print("出现警告")
# 		print(w)
# 	}, error = function(e){
# 		print("出现错误")
# 		print(e)
# 	},finally = {
# 		# 这是运行正常时，应该怎么做，可以用print打印出来，也可以执行其它命令
# 		print("good")
# 		pdf(file.path("C:/Users/xi.xu/Desktop/热图决策树--小工具", 'test15.pdf'), width = 12, height = 6, onefile=FALSE)
# 		heat_tree(xx_df, target_lab = 'sex',show = 'heat-tree')
# 		dev.off()
# 	}
# )
