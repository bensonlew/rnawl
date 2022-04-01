#!/usr/bin/env Rscript
times<-Sys.time()
library('getopt');
options(bitmapType='cairo')
options(scipen = 200)
spec = matrix(c(
	'infile','i',0,'character',
	'outfile','o',0,'character',
	'winsize','w',0,'character',
	'method','m',0,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--infile	the input hapmap file
	--outfile	the trait file 
	--winsize	the window size
	--method	the sliding method default bp or num
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) )   { print_usage(spec) }
if ( is.null(opt$infile))   { print_usage(spec)}
if ( is.null(opt$outfile))  { print_usage(spec) }
if ( is.null(opt$winsize))  {opt$winsize=1000000}
if (is.null(opt$method)){opt$method="bp"}
opt$winsize=as.numeric(opt$winsize)
library(dplyr)
library("Rcpp")
argv <- commandArgs(trailingOnly=F)
temp <- '--file='
script.abspath <- sub(temp, "", argv[grep(temp, argv)])
script.dirname <- dirname(script.abspath)
sourceCpp(paste(script.dirname,"countSNPs.cpp",sep="/"))
tricubeStat <- function(POS, Stat, windowSize = 2e6, ...)
{
    if (windowSize <= 0)
        stop("A positive smoothing window is required")
    stats::predict(locfit::locfit(Stat ~ locfit::lp(POS, h = windowSize, deg = 0),maxk=100000), POS)
}
getG <- function(LowRef, HighRef, LowAlt, HighAlt)
{
    exp <- c(
        (LowRef + HighRef) * (LowRef + LowAlt) / (LowRef + HighRef + LowAlt + HighAlt),
        (LowRef + HighRef) * (HighRef + HighAlt) / (LowRef + HighRef + LowAlt + HighAlt),
        (LowRef + LowAlt) * (LowAlt + HighAlt) / (LowRef + HighRef + LowAlt + HighAlt),
        (LowAlt + HighAlt) * (HighRef + HighAlt) / (LowRef + HighRef + LowAlt + HighAlt)
    )
    obs <- c(LowRef, HighRef, LowAlt, HighAlt)

    G <-
        2 * (rowSums(obs * log(
            matrix(obs, ncol = 4) / matrix(exp, ncol = 4)
        )))
    return(G)
}
psites<-function(d){
	c(1:length(d))
}

data<-read.table(opt$infile,head=TRUE,comment.char="^")
data$depth=data$n1+data$n2+data$n3+data$n4
if(opt$method == "num"){
	data<-data %>%
		group_by(X.chr) %>%
		mutate(psites = psites(pos) )
}else{
	data$psites=data$pos
}
data$G=getG(LowRef = data$n2,HighRef = data$n1,LowAlt = data$n4,HighAlt = data$n3)
data$G[data$G == "NaN"]=0.00000000001
print("haha")
data <- data %>%
	group_by(X.chr) %>%
	mutate(Gprime = tricubeStat(POS = psites, Stat = G, opt$winsize))
data<-data %>% group_by(X.chr) %>%
	mutate(ncount = countSNPs_cpp(POS = psites, windowSize = opt$winsize))
data$Gprime[data$ncount < 10]=mean(data$Gprime)
lnGprime <- log(data$Gprime)
medianLogGprime <- median(lnGprime)
MAD <-median(medianLogGprime - lnGprime[lnGprime <= medianLogGprime])
trimGprime <-data$Gprime[lnGprime - median(lnGprime) <= 5.2 * MAD]
medianTrimGprime <- median(trimGprime)
modeTrimGprime <- modeest::mlv(x = trimGprime, bw = 0.5, method = "hsm")
muE <- log(medianTrimGprime)
varE <- abs(muE - log(modeTrimGprime))
data$Gpval <-1 - plnorm(q = data$Gprime,meanlog = muE,sdlog = sqrt(varE))
data$qvalue = p.adjust(p = data$Gpval, method = "BH")
data=arrange(data,X.chr,pos)
write.table(file=opt$outfile,data,quote=F,row.names=F)

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
q()
