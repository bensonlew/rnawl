# # # # ~/app/bioinfo/rna/miniconda2/bin/R
# # export PATH=~/app/library/curl-7.54/bin/:$PATH
# # export PATH=/mnt/lustre/users/sanger-dev/app/gcc/7.2.0/bin/:$PATH
# # export LD_LIBRARY_PATH=/mnt/lustre/users/sanger-dev/app/gcc/7.2.0/lib64:$LD_LIBRARY_PATH

# options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# install.packages("https://cran.r-project.org/src/contrib/Archive/bitops/bitops_1.0-6.tar.gz", repos=NULL, type="source")
# install.packages("https://cran.r-project.org/src/contrib/Archive/caTools/caTools_1.17.1.2.tar.gz", repos=NULL, type="source")
# install.packages("gtools")
# install.packages("https://mirrors.tuna.tsinghua.edu.cn/CRAN/src/contrib/gplots_3.1.1.tar.gz", repos=NULL, type="source")
# install.packages("flexclust")
# install.packages("additivityTests")
# ######## install.packages("https://mirrors.tuna.tsinghua.edu.cn/CRAN/src/contrib/biclust_2.0.3.tar.gz", repos=NULL, type="source")
# install.packages("https://cran.r-project.org/src/contrib/Archive/biclust/biclust_2.0.2.tar.gz", repos=NULL, type="source")
# install.packages("https://mirrors.tuna.tsinghua.edu.cn/CRAN/src/contrib/proxy_0.4-26.tar.gz", repos=NULL, type="source")
# install.packages("https://mirrors.tuna.tsinghua.edu.cn/CRAN/src/contrib/e1071_1.7-7.tar.gz", repos=NULL, type="source")
# install.packages("https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.1.tar.gz", repos=NULL, type="source")
# install.packages("https://cran.r-project.org/src/contrib/Archive/latticeExtra/latticeExtra_0.6-28.tar.gz", repos=NULL, type="source")
# install.packages("testthat")
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# install.packages("RCurl")
# BiocManager::install("GenomeInfoDb")
# BiocManager::install("GenomicRanges")
# BiocManager::install("SummarizedExperiment")
# BiocManager::install("Hmisc")

# BiocManager::install("edgeR")
# BiocManager::install("QUBIC")
# BiocManager::install("runibic")

# install.packages("XML")
# BiocManager::install("DESeq2")


library('getopt');

spec = matrix(c(
	'infile','i',1,'character',
	'outdir','o',1,'character',
	'method','m',2,'character',
	'preprocess','p',2,'logical',
	'heatmap','d',2,'logical',
	'ngenes','n',2,'integer',
	'missing','s',2,'character',
	'format','f',2,'integer',
	'transform','t',2,'integer',
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
	--method(-m)	 the method to biocluster(default BCCC()) : [BCCC(),BCQU(),BCUnibic(),BCXmotifs(),BCPlaid(),BCSpectral(),BCBimax(),BCQuest()]
	--preprocess(-p)	Will Do Pre-process before bicluster.(default FALSE)
	--heatmap(-d)	Will Draw heatmap of all subcluster.(default FALSE)
	--ngenes(-n)	Most variable genes to include. (default 1000)
	--missing(-s)	Missing values imputation. (default geneMedian)
	--format(-f)	Input_data File Format. (default 1,reads_count)
	--transform(-t)	Transform counts data for clustering. (default 1,EdgeR log2(CPM+c))
	--help(-h)		usage
\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$infile)) { print_usage(spec)}
if ( is.null(opt$outdir)) { print_usage(spec)}
if ( is.null(opt$method)){ opt$method="BCCC()"}
if ( is.null(opt$preprocess)){ opt$preprocess=FALSE}else{opt$preprocess=TRUE}
if ( is.null(opt$heatmap)){ opt$heatmap=FALSE }else{opt$heatmap=TRUE}
if ( is.null(opt$ngenes)){ opt$ngenes=1000}
if ( is.null(opt$missing)){ opt$missing="geneMedian"}
if ( is.null(opt$format)){ opt$format=1}
if ( is.null(opt$transform)){ opt$transform=1}

inFile = opt$infile
output_dir = opt$outdir
input_biclustMethod = opt$method
preprocess = opt$preprocess
paint_heatmap = opt$heatmap
input_nGenesBiclust = opt$ngenes
input_missingValue = opt$missing
input_dataFileFormat = opt$format  #2,3(FPKM,other program) not test yet 
input_CountsTransform = opt$transform #2,3(VST,rlog) not test yet 

print(paste("inFile： ",inFile))
print(paste("output_dir： ",output_dir))
print(paste("input_biclustMethod： ",input_biclustMethod))
print(paste("preprocess： ",preprocess))
print(paste("paint_heatmap： ",paint_heatmap))
print(paste("input_nGenesBiclust： ",input_nGenesBiclust))
print(paste("input_missingValue： ",input_missingValue))
print(paste("input_dataFileFormat： ",input_dataFileFormat))
print(paste("input_CountsTransform： ",input_CountsTransform))


library("e1071",verbose=FALSE)
library(biclust,verbose=FALSE)


detectGroups <- function (x){  # x are col names
  tem <- gsub("[0-9]*$","",x) # Remove all numbers from end
  #tem = gsub("_Rep|_rep|_REP","",tem)
  tem <- gsub("_$","",tem); # remove "_" from end
  tem <- gsub("_Rep$","",tem); # remove "_Rep" from end
  tem <- gsub("_rep$","",tem); # remove "_rep" from end
  tem <- gsub("_REP$","",tem)  # remove "_REP" from end
  return( tem )
}


readData <- function(inFile ) {
  # these packages moved here to reduce loading time
  library(edgeR,verbose=FALSE) # count data D.E.
  library(DESeq2,verbose=FALSE) # count data analysis
  dataTypeWarning =0
  dataType =c(TRUE)
  #---------------Read file
  x <- read.csv(inFile)	# try CSV
  if(dim(x)[2] <= 2 )   # if less than 3 columns, try tab-deliminated
    x <- read.table(inFile, sep="\t",header=TRUE)
  #-------Remove non-numeric columns, except the first column
  for(i in 2:dim(x)[2])
    dataType = c( dataType, is.numeric(x[,i]) )
  if(sum(dataType) <=2) return (NULL)  # only less than 2 columns are numbers
  x <- x[,dataType]  # only keep numeric columns
  # rows with all missing values
  ix = which( apply(x[,-1],1, function(y) sum( is.na(y) ) ) != dim(x)[2]-1 )
  x <- x[ix,]
  dataSizeOriginal = dim(x); dataSizeOriginal[2] = dataSizeOriginal[2] -1
  x[,1] <- toupper(x[,1])
  x[,1] <- gsub(" ","",x[,1]) # remove spaces in gene ids
  x = x[order(- apply(x[,2:dim(x)[2]],1,sd) ),]  # sort by SD
  x <- x[!duplicated(x[,1]) ,]  # remove duplicated genes
  rownames(x) <- x[,1]
  x <- as.matrix(x[,c(-1)])
  # remove "-" or "." from sample names
  colnames(x) = gsub("-","",colnames(x))
  colnames(x) = gsub("\\.","",colnames(x))
  #cat("\nhere",dim(x))
  # missng value for median value
  if(sum(is.na(x))>0) {# if there is missing values
    if(input_missingValue =="geneMedian") {
      rowMedians <- apply(x,1, function (y)  median(y,na.rm=T))
      for( i in 1:dim(x)[2] ) {
        ix = which(is.na(x[,i]) )
        x[ix,i] <- rowMedians[ix]
      }
    } else if(input_missingValue =="treatAsZero") {
      x[is.na(x) ] <- 0
    } else if (input_missingValue =="geneMedianInGroup") {
      sampleGroups = detectGroups( colnames(x))
      for (group in unique( sampleGroups) ){
        samples = which( sampleGroups == group )
        rowMedians <- apply(x[,samples, drop=F],1, function (y)  median(y,na.rm=T))
        for( i in  samples ) {
          ix = which(is.na(x[ ,i] ) )
          if(length(ix) >0 )
            x[ix, i  ]  <- rowMedians[ix]
        }
      }
      # missing for entire sample group, use median for all samples
      if(sum(is.na(x) )>0 ) {
        rowMedians <- apply(x,1, function (y)  median(y,na.rm=T))
        for( i in 1:dim(x)[2] ) {
          ix = which(is.na(x[,i]) )
          x[ix,i] <- rowMedians[ix]
        }
      }
    }
  }
  # Compute kurtosis
  mean.kurtosis = mean(apply(x,2, kurtosis),na.rm=T)
  rawCounts = NULL
  pvals= NULL
  if (input_dataFileFormat == 2 ) {  # if FPKM, microarray
    if ( is.integer(x) ) dataTypeWarning = 1;  # Data appears to be read counts
    #-------------filtering
    #tem <- apply(x,1,max)
    #x <- x[which(tem > input_lowFilter),]  # max by row is at least
    x <- x[ which( apply( x, 1,  function(y) sum(y >= input_lowFilter)) >= input_NminSamples2 ) , ]
    x <- x[which(apply(x,1, function(y) max(y)- min(y) ) > 0  ),]  # remove rows with all the same levels
    #--------------Log transform
    # Takes log if log is selected OR kurtosis is big than 100
    if ( (input_transform == TRUE) | (mean.kurtosis > kurtosis.log ) )
      x = log(x+abs( input_logStart),2)
    tem <- apply(x,1,sd)
    x <- x[order(-tem),]  # sort by SD
  } else
    if( input_dataFileFormat == 1) {  # counts data
      # data not seems to be read counts
      if(!is.integer(x) & mean.kurtosis < kurtosis.log ) {
        dataTypeWarning = -1
      }
      # not used as some counts data like those from CRISPR screen
      #validate(   # if Kurtosis is less than a threshold, it is not read-count
      #	need(mean.kurtosis > kurtosis.log, "Data does not seem to be read count based on distribution. Please double check.")
      # )
      x <- round(x,0) # enforce the read counts to be integers. Sailfish outputs has decimal points.
      #x <- x[ which( apply(x,1,max) >= input_minCounts ) , ] # remove all zero counts
      # remove genes if it does not at least have minCounts in at least NminSamples
      #x <- x[ which( apply(x,1,function(y) sum(y>=input_minCounts)) >= input_NminSamples ) , ]  # filtering on raw counts
      # using counts per million (CPM) for filtering out genes.
      # CPM matrix                  #N samples > minCounts
      x <- x[ which( apply( cpm(DGEList(counts = x)), 1,
                            function(y) sum(y>=input_minCounts)) >= input_NminSamples ) , ] # 0.5 ， 1
      rawCounts = x; # ???
      # construct DESeqExpression Object
      # colData = cbind(colnames(x), as.character(detectGroups( colnames(x) )) )
      tem = rep("A",dim(x)[2]); tem[1] <- "B"   # making a fake design matrix to allow process, even when there is no replicates
      colData = cbind(colnames(x), tem )
      colnames(colData)  = c("sample", "groups")
      dds <- DESeqDataSetFromMatrix(countData = x, colData = colData, design = ~ groups)
      dds <- estimateSizeFactors(dds) # estimate size factor for use in normalization later for started log method
      # regularized log  or VST transformation
      if( input_CountsTransform == 3 ) { # rlog is slow, only do it with 10 samples
        if(dim(x)[2]<=20 ) {
          x <- rlog(dds, blind=TRUE); x <- assay(x) } else
            x <- log2( counts(dds, normalized=TRUE) + input_countsLogStart ) # 4 default
      }
      else {
        if ( input_CountsTransform == 2 ) {    # vst is fast but aggressive
          x <- vst(dds, blind=TRUE)
          x <- assay(x)
        } else{  # normalized by library sizes and add a constant.
          x <- log2( counts(dds, normalized=TRUE) + input_countsLogStart )   # log(x+c)   # 4 default
          # This is equivalent to below. But the prior counts is more important
          #x <- cpm(DGEList(counts = x),log=TRUE, prior.count=input_countsLogStart )  #log CPM from edgeR
          #x <- x-min(x)  # shift values to avoid negative numbers
        }
      }
    } else
      if( input_dataFileFormat == 3)	{  # other data type
        #neg_lfc neg_fdr pos_lfc pos_fdr
        #11       1      11       1
        n2 = ( dim(x)[2] %/% 2) # 5 --> 2
        # It looks like it contains P values
        # ranges of columns add 0.2 and round to whole. For P value columns this should be 1
        tem = round( apply(x, 2, function( y) max(y)- min(y))  + .2)
        if( sum(tem[(1:n2)*2  ] ==  1 ) == n2 |
            sum(tem[(1:n2)*2-1  ] ==  1 ) == n2 ) {
          x = x[,1:(2*n2) ,drop=FALSE ] # if 5, change it to 4
          if(tem[2] == 1) { # FDR follows Fold-change
            pvals = x [,2*(1:n2 ),drop=FALSE ]  # 2, 4, 6
            x = x[, 2*(1:n2 )-1,drop=FALSE]   # 1, 3, 5
          } else {	# FDR follows Fold-change
            pvals = x [,2*(1:n2 )-1,drop=FALSE ]  # 2, 4, 6
            x = x[, 2*(1: n2 ),drop=FALSE]   # 1, 3, 5
          }
        }
        ix =  which(apply(x,1, function(y) max(y)- min(y) ) > 0  )
        x <- x[ix,]  # remove rows with all the same levels
        if(!is.null(pvals) )
          pvals = pvals[ix,]
      }
  dataSize = dim(x);
  if(!(dim(x)[1]>5 & dim(x)[2]>1))
    stop ( "Data file not recognized. Please double check.")
  sampleInfoDemo=NULL
  if( input_goButton >0)
    sampleInfoDemo <- t( read.csv(demoDataFile2,row.names=1,header=T,colClasses="character") )
  finalResult <- list(data = as.matrix(x))
  return(finalResult)
}


biclustering <- function( ){
  library(biclust,verbose=FALSE)
  if(input_biclustMethod == "BCQU()" )
    library(QUBIC,verbose=FALSE) # have trouble installing on Linux
  if(input_biclustMethod == "BCUnibic()" )
    library(runibic,verbose=FALSE) #
  x <- convertedData.out
  n=input_nGenesBiclust
  if(n>dim(x)[1]) n = dim(x)[1] # max	as data
  if(n<10) n = 10
  if(n> 50000 ) n = 50000
  x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
  if( input_biclustMethod == "BCXmotifs()"  )
    x<-discretize(x)
  #res <- biclust::biclust(as.matrix( x), method = BCQU())
  runR = paste( "res <- biclust::biclust(as.matrix( x), method =", input_biclustMethod ,")" )
  eval(parse(text = runR ) )
  return(list( x=x, res=res)	)
}


biclustHeatmap <- function (input_selectBicluster){
  res = biclustering()$res
  if( res@Number == 0 ) { plot.new(); text(0.5,0.5, "No cluster found!")} else {
    
    x = biclust::bicluster(biclustering.out$x, res, as.numeric( input_selectBicluster)  )[[1]]
    
    par(mar = c(5, 4, 1.4, 0.2))
    
    if( dim(x)[1] <=30 )
      heatmap.2(x,  Rowv =T,Colv=F, dendrogram ="none",
                col=heatColors[as.integer(input_heatColors1),], density.info="none", trace="none", scale="none", keysize=.2
                ,key=F, #labRow = T,
                ,margins = c(8, 24)
                ,cexRow=1
                #,srtCol=45
                ,cexCol=1.  # size of font for sample names
      ) else 		
        heatmap.2(x,  Rowv =T,Colv=F, dendrogram ="none",
                  col=heatColors[as.integer(input_heatColors1),], density.info="none", trace="none", scale="none", keysize=.2
                  ,key=F, labRow = F,
                  ,margins = c(8, 4)
                  ,cexRow=1
                  #,srtCol=45
                  ,cexCol=1.  # size of font for sample names
        )				
  }
}




input_minCounts = 0.5
input_NminSamples = 1
input_lowFilter = -1000
input_NminSamples2 = 1
input_transform = FALSE
input_logStart = 1
kurtosis.log = 50
input_countsLogStart = 4
input_goButton = 0 #do not change it
input_heatColors1 = 2
# input_dataFileFormat = 1  #2,3(FPKM,other program) not test yet 
# input_CountsTransform = 1 #2,3(VST,rlog) not test yet 
# input_missingValue ="geneMedian"
# input_nGenesBiclust = 2000
# input_biclustMethod = "BCCC()"


# inFile = "C:/Users/xi.xu/Desktop/双聚类小工具/counts_cmp_use_small_good.txt"
# output_dir = "C:/Users/xi.xu/Desktop/双聚类小工具/"
# preprocess = TRUE
# paint_heatmap = TRUE


if (preprocess == TRUE) {
	library(edgeR,verbose=FALSE) # count data D.E.
	library(DESeq2,verbose=FALSE) # count data analysis
	readData.out = readData(inFile = inFile)
	convertedData.out = readData.out$data
}else{
	convertedData.out = as.matrix(read.table(inFile,header=T,row.names= 1))
}

biclustering.out = biclustering()
write.table(biclustering.out$res@Number, file.path(output_dir, "bicluster.num"),row.names=FALSE, col.names=FALSE)
if( biclustering.out$res@Number > 0 ) {
	for ( i in 1:biclustering.out$res@Number) {
		biclustering.out_res_1 = bicluster(biclustering.out$x, biclustering.out$res, i )[[1]]
		write.table(biclustering.out_res_1, file=file.path(output_dir, paste("sub_cluster_",as.character(i),".tsv",sep="")), sep="\t",row.names=TRUE, col.names=NA, quote = FALSE)
	}
}

if (paint_heatmap == TRUE){
	if( biclustering.out$res@Number > 0 ) {
		library(gplots,verbose=FALSE)
		#
		hmcols <- colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","#E0F3F8", "#91BFDB", "#4575B4")))(75)
		heatColors = rbind(greenred(75), bluered(75), colorpanel(75,"green","black","magenta"),colorpanel(75,"blue","yellow","red"),hmcols)
		rownames(heatColors) = c("Green-Black-Red","Blue-White-Red","Green-Black-Magenta","Blue-Yellow-Red","Blue-white-brown")
		colorChoices = setNames(1:dim(heatColors)[1],rownames(heatColors)) # for pull down menu
		#
		for ( i in 1:biclustering.out$res@Number) {
			pdf(file = file.path(output_dir, paste("sub_cluster_",as.character(i),".pdf",sep="")))
			biclustHeatmap(input_selectBicluster=i)
			dev.off()
			# svg(file = file.path(output_dir, paste("sub_cluster_",as.character(i),".svg",sep="")))
			# biclustHeatmap(input_selectBicluster=i)
			# dev.off()
		}
	}
}


# input_biclustMethod = "BCXmotifs()"
# biclustering.out = biclustering()
# biclustering.out$res

# input_biclustMethod = "BCUnibic()"
# biclustering.out = biclustering()
# biclustering.out$res

# input_biclustMethod = "BCQU()"
# biclustering.out = biclustering()
# biclustering.out$res

# input_biclustMethod = "BCPlaid()"
# biclustering.out = biclustering()
# biclustering.out$res

# input_biclustMethod = "BCSpectral()"
# biclustering.out = biclustering()
# biclustering.out$res

# input_biclustMethod = "BCBimax()"
# biclustering.out = biclustering()
# biclustering.out$res

# input_biclustMethod = "BCQuest()"
# biclustering.out = biclustering()
# biclustering.out$res

# input_biclustMethod = "BCCC()"
# biclustering.out = biclustering()
# biclustering.out$res