times<-Sys.time()
library('getopt');
library('Cairo')
options(bitmapType='cairo')
spec = matrix(c(
	'mark','m',1,'character',
	'trt','t',1,'character',
	'out','o',1,'character',
	'num','n',1,'character',
	'pop','p',1,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4)
opt = getopt(spec)
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript Rqtl-NOCP.r --mark  --out --num --pop
	
Usage:
	--mark	map file
	--out	out dir
	--help		usage
\n")
	q(status=1);
}
times<-Sys.time()
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$mark) ) { print_usage(spec) }
if ( is.null(opt$out) ) { print_usage(spec) }

d<-read.table(opt$mark,head=TRUE,sep=",")
colnames(d)[2:3]=c("chr","pos")
chr=unique(d$chr);
for(i in chr){
	pos=(1:length(d$pos[d$chr==i]))
	d$start[d$chr==i]=pos
	d$end[d$chr==i]=pos+1
}
print(unique(d$chr))

cols=c("red","blue","orange","grey40")
names(cols)=c("A","B","H","-")
chr.len <- tapply(d$end, d$chr, max)
x.brks <- cumsum(chr.len) - chr.len/2
chr.cum.len <- c(0, cumsum(chr.len)[-length(chr.len)])
d$start<-d$start+chr.cum.len[d$chr]
d$end<-d$end+chr.cum.len[d$chr]
names(chr.cum.len) <- names(chr.len)

# pdf(paste(opt$out,"pdf",sep="."),height=900,width=1600);
# plot(5, type = "n", xlim = c(0, max(d$end,na.rm=TRUE)), ylim = c(0,ncol(d) - 4), xaxt = "n",main="Haplotype")
# axis(side = 1, at = x.brks, labels = unique(d$chr))
# for (i in 4:ncol(d)) {
	# if(colnames(d)[i]=="start"|| colnames(d)[i]=="end"){next;}
    # rect(xleft = d$start, ybottom = i - 5, xright = d$end,ytop = i - 4.2, col = cols[as.character(d[,i])], border = NA)
	# rect(xleft = d$start, ybottom = i - 4.2, xright = d$end,ytop = i - 4, col = cols[rep("-",length(as.character(d[,i])))], border = NA)

# }
# abline(v = c(0, cumsum(chr.len)), col = "grey80", lwd = 0.5)
# dev.off()

CairoPNG(paste(opt$out,"png",sep="."),height=900,width=1600);
plot(5, type = "n", xlim = c(0, max(d$end,na.rm=TRUE)), ylim = c(0,ncol(d) - 4), xaxt = "n",main="Haplotype")
axis(side = 1, at = x.brks, labels = unique(d$chr))
for (i in 4:ncol(d)) {
	if(colnames(d)[i]=="start"|| colnames(d)[i]=="end"){next;}
    rect(xleft = d$start, ybottom = i - 5, xright = d$end,ytop = i - 4.2, col = cols[as.character(d[,i])], border = NA)
	rect(xleft = d$start, ybottom = i - 4.2, xright = d$end,ytop = i - 4, col = cols[rep("-",length(as.character(d[,i])))], border = NA)
}
abline(v = c(0, cumsum(chr.len)), col = "grey80", lwd = 0.5)
dev.off()
for (i in chr){
	subd<-d[d$chr==i,];
	pos=c(1:length(subd$pos[subd$chr==i]))
	subd$start[subd$chr==i]=pos
	subd$end[subd$chr==i]=pos+1
	# pdf(paste(opt$out,i,"pdf",sep="."),height=900,width=1600);
	# plot(5, type = "n", xlim = c(0, max(subd$end,na.rm=TRUE)), ylim = c(0,ncol(subd) - 4),main="Haplotype",xlab=paste("LG",i,sep=" "),ylab="Samples")
	# for (j in 4:ncol(subd)) {
		# if(colnames(subd)[j]=="start"|| colnames(subd)[j]=="end"){next}
		# rect(xleft = subd$start, ybottom = j - 5, xright = subd$end,ytop = j - 4.2, col = cols[as.character(subd[,j])], border = NA)
		# rect(xleft = subd$start, ybottom = j - 4.2, xright = subd$end,ytop = j - 4, col = cols[rep("-",length(as.character(subd[,j])))], border = NA)
	# }
	# dev.off()
	CairoPNG(paste(opt$out,i,"png",sep="."),height=900,width=1600);
	plot(5, type = "n", xlim = c(0, max(subd$end,na.rm=TRUE)), ylim = c(0,ncol(subd) - 4),main="Haplotype",xlab=paste("LG",i,sep=" "),ylab="Samples")
	for (j in 4:ncol(subd)) {
		if(colnames(subd)[j]=="start"|| colnames(subd)[j]=="end"){next}
		rect(xleft = subd$start, ybottom = j - 5, xright = subd$end,ytop = j - 4.2, col = cols[as.character(subd[,j])], border = NA)
		rect(xleft = subd$start, ybottom = j - 4.2, xright = subd$end,ytop = j - 4, col = cols[rep("-",length(as.character(subd[,j])))], border = NA)
	}
	dev.off()

}

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
