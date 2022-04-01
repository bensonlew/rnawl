#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
my %opts;
my $VERSION = "V1_20160510";

GetOptions( \%opts,"i=s","o=s","g=s","s=s","group=s","t=s","w=i","h=i","cstar=i","cel=i","axecel=s","clabel=i","cpoint=i","grid=s","addaxes=s","csub=i","possub=s","cgrid=i","l=s", "d=s");

my $usage = <<"USAGE";
       Program : $0   
       Discription:   
       Version : $VERSION
       Contact : hao.gao\@majorbio.com
       Last_modified : zouxuan 2017.11.09
       Usage :perl $0 [options]		
			-i	    * input genus table file 
                        -g          * clusters of input eg:1,2,3,...
                        -s          * One of the best species relative abundance of the cluster eg:\\"Prevotella\\",\\"Enterococcus,Escherichia-Shigella,Prevotella\\",\\"Bacteroides\\",....
		        -o	    * output dir
                        -group      * group of samples eg: #name  group
#                        -t          *  type of plot eg:PCOA or BCA    由于BCA与pcoa结果都生成，因此暂取消该参数  zouxuan
                        -d        dist_matrix #add by zouxuan 20171120
                        -w	      width of pdf  default :10
                        -h            height of pdf  default:
                        -cstar        a number between 0 and 1 which defines the length of the star size eg: [0,1] eg��default��1
                        -cel               the ellipse size eg��default��1.5
                        -axecel            a logical value indicating whether the ellipse axes should be drawn  eg��default��T
                        -clabel            a character size for star  eg��default��1
                        -cpoint            a character size for plotting the points eg��default��1
                        -grid              the background of the plot should be drawn  eg��default��T
                        -addaxes           the axes should be plotted  eg��default��T
                        -csub              a character size for the legend  eg��default��1
                        -possub            a string of characters indicating the sub-title position ("topleft", "topright", "bottomleft", "bottomright") eg��default��"bottomleft"
                        -cgrid             the mesh of the grid
                        -l                 whether the sample name should be draw  eg:default:F
       Example:perl ../plot-Enterotyping.pl -i otu_3600.txt -g 1,2,3 -s \"Prevotella\",\"Enterococcus,Escherichia-Shigella,Prevotella\",\"Bacteroides\" -o gao -group ../group -t BCA -l T 

USAGE

die $usage if(!($opts{i}&&$opts{o}&&$opts{s}&&$opts{g}&&$opts{group}));

$opts{w}=defined $opts{w}?$opts{w}:10;
$opts{h}=defined $opts{h}?$opts{h}:10;
$opts{cstar}=defined $opts{cstar}?$opts{cstar}:1;
$opts{cel}=defined $opts{cel}?$opts{cel}:1.5;
$opts{axecel}=defined $opts{axecel}?$opts{axecel}:"T";
$opts{clabel}=defined $opts{clabel}?$opts{clabel}:1;
$opts{cpoint}=defined $opts{cpoint}?$opts{cpoint}:1;
$opts{grid}=defined $opts{grid}?$opts{grid}:"F";
$opts{addaxes}=defined $opts{addaxes}?$opts{addaxes}:"T";
$opts{csub}=defined $opts{csub}?$opts{csub}:1;
$opts{possub}=defined $opts{possub}?$opts{possub}:"bottomleft";
$opts{cgrid}=defined $opts{cgrid}?$opts{cgrid}:1;
$opts{l}=defined $opts{l}?$opts{l}:"F";
$opts{d}=defined $opts{d}?$opts{d}:"none"; #add by zouxuan 20171120
### plot ###

open RCMD, ">plot-Enterotyping.r";
print RCMD "

library(\"cluster\")
library(\"clusterSim\")
library(\"ade4\")
data=read.delim(\"$opts{i}\",sep=\"\\t\",header=T, row.names=1, check.names=F)
# data=data[-1,]
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
 
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
 
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) {
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, \"method\") <- \"dist\"
  return(resultsMatrix)
}
if(\"$opts{d}\"  == \"none\"){
    data.dist=dist.JSD(data)
}else{
   file = read.table(\"$opts{d}\",sep=\"\\t\",header=T, row.names=1, check.names=F)
   data.dist = as.dist(file)
}#add by zouxuan 20171120

pam.clustering=function(x,k) { 
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)\$clustering)
  return(cluster)
}
nclusters = NULL
for (k in 1:9) {
  if (k==1) {
    nclusters[k]=NA
  } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k]=index.G1(t(data),data.cluster_temp, d = data.dist, centrotypes = \"medoids\")
  }
}
max=which.max(nclusters)
data.cluster=pam.clustering(data.dist, k=max)

###��ͼ###
mycol <-c(\"#CD0000\",\"#3A89CC\",\"#769C30\",\"#D99536\",\"#7B0078\",\"#BFBC3B\",\"#6E8B3D\",\"#00688B\",\"#C10077\",\"#CAAA76\",\"#EEEE00\",\"#458B00\",\"#8B4513\",\"#008B8B\",\"#6E8B3D\",\"#8B7D6B\",\"#7FFF00\",\"#CDBA96\",\"#ADFF2F\")
mypch <-c(21:25,3,4,7,9,8,10,15:18,0:14)
col<-mycol[1:max]
map<-read.table(\"$opts{group}\",sep=\"\",check.names=F)
pch<-mypch[1:length(unique(map[,2]))]
group<-unique(map[,2])
nem<-names(data)

for (i in 1:length(unique(map[,2]))){
map[which(map[,2]==group[i]),3]=pch[i]
}
id<-which(map[,1] \%in\% nem)
pchl<-map[which(map[,1] \%in\% nem),3]
name=names(data[0,])
#name<-factor(name)
group1=data.cluster
all=cbind(name,group1)
map\$v4=all[,2]
zp<-max(map\$v4)
zp2<-sum(all[,1]==map[,1])
zp2<-sum(which(all[,1]==map[,1]))
zhangpeng1<-length(map[,1])
map\$v5<-rep(1,zhangpeng1)
zhangpeng2<-tapply(map\$v5,list(map\$v4,map\$V2),sum)
zhangpeng3<-colnames(zhangpeng2)
zhangpeng4<-length(zhangpeng2[,1])
zhangpeng6<-length(zhangpeng2[1,])
zhangpeng11<-c(1:zhangpeng4)
zhangpeng12<-zhangpeng2[,1:zhangpeng6]
zhangpeng5<-cbind(zhangpeng11,zhangpeng12)
colnames(zhangpeng5)<-c(\"Enterotype\",zhangpeng3)
zhangpeng5[is.na(zhangpeng5)] <- 0

file_1 <- paste(\"$opts{o}/\",\"summary.txt\",sep=\"\")
write.table(zhangpeng5,file=file_1,sep=\"\t\",quote = F,row.names=F,col.names=T)

result <- vector()
groups <- unique(map[,2])
obs.pca=dudi.pca(data.frame(t(data)), scannf=F, nf=10)
obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=2)
obs.pcoa=dudi.pco(data.dist, scannf=F, nf=2)
a<-c($opts{s})
names(a)<-c($opts{g})
b<-a[data.cluster]
names(b) <- NULL


fac2disj<- function(fac, drop = FALSE) {
  ## Returns the disjunctive table corrseponding to a factor
  n <- length(fac)
  fac <- as.factor(fac)
  if(drop)
    fac <- factor(fac)
  x <- matrix(0, n, nlevels(fac))
  x[(1:n) + n * (unclass(fac) - 1)] <- 1
  dimnames(x) <- list(names(fac), as.character(unique(fac)))
  return(data.frame(x, check.names = FALSE))
}

#if(\"$opts{l}\" == \"T\"){
#s.class(obs.bet\$ls, fac=as.factor(b), pch=pchl,col=col,sub=\"Between-class analysis\",cstar=$opts{cstar},cellipse=$opts{cel},axesell=$opts{axesel},clabel=$opts{clabel},cpoint=$opts{cpoint},grid=$opts{grid},addaxes=$opts{addaxes},csub=$opts{csub},possub=\"$opts{possub}\",cgrid=$opts{cgrid})
#########zhangpeng########
xax = 1
yax = 2
pch=pchl
col=col
sub=\"Between-class analysis\"
cstar=$opts{cstar}
cellipse=$opts{cel}
axesell=$opts{axesel}
clabel=$opts{clabel}
cpoint=$opts{cpoint}
grid=$opts{grid}
addaxes=$opts{addaxes}
csub=$opts{csub}
possub=\"$opts{possub}\"
cgrid=$opts{cgrid}
xlim = NULL
ylim = NULL
origin = c(0, 0)
include.origin = TRUE
pixmap = NULL
contour = NULL
area = NULL
add.plot = FALSE


for (select in c(\"BCA\",\"pcoa\")) {
    if(select == \"BCA\"){
        my_colum = obs.bet\$ls
    }else{
    my_colum = obs.pcoa\$li}
    pdf_file = paste(select,\".pdf\",sep = \"\")
    pdf(pdf_file,$opts{w},$opts{h})
    opar <- par(mar = par(\"mar\"))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    on.exit(par(opar))
    file2 <- paste(\"$opts{o}/\",paste(select,\"_point.txt\",sep=\"\"),sep=\"\")
    if(ncol(my_colum)==1){
    col_name = c(\"Sample_ID\",\"x\")
    write.table(cbind(rownames(my_colum),my_colum),file=file2,sep=\"\t\",quote = F,row.names=F,col.names=col_name)
    next
    }else{
    dfxy <- data.frame(my_colum)
    write.table(cbind(rownames(my_colum),my_colum),file=file2,sep=\"\t\",quote = F,row.names=F,col.names=c(\"Sample_ID\",\"x\",\"y\"))
    }
    if (!is.data.frame(my_colum))
        stop(\"Non convenient selection for dfxy\")
    if (any(is.na(my_colum)))
        stop(\"NA non implemented\")
    if (!is.factor(as.factor(b)))
        stop(\"factor expected for fac\")
    dfdistri <- fac2disj(as.factor(b)) * rep(1, length(as.factor(b)))
    coul <- col
    w1 <- unlist(lapply(dfdistri, sum))
    dfdistri <- t(t(dfdistri)/w1)
    coox <- as.matrix(t(dfdistri)) \%*\% dfxy[, xax]
    cooy <- as.matrix(t(dfdistri)) \%*\% dfxy[, yax]
    if (nrow(dfxy) != nrow(dfdistri))
        stop(paste(\"Non equal row numbers\", nrow(dfxy), nrow(dfdistri)))
    coo <- scatterutil.base(dfxy = dfxy, xax = xax, yax = yax,
        xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes,
        cgrid = cgrid, include.origin = include.origin, origin = origin,
        sub = sub, csub = csub, possub = possub, pixmap = pixmap,
        contour = contour, area = area, add.plot = add.plot)
    if (cpoint > 0)
        for (i in 1:ncol(dfdistri)) {
            pch <- rep(pch, length = nrow(dfxy))
            points(coo\$x[dfdistri[, i] > 0], coo\$y[dfdistri[,
                i] > 0], pch = pch[dfdistri[, i] > 0], cex = par(\"cex\") *
                cpoint, col = coul[i])
        }
    if (cstar > 0)
        for (i in 1:ncol(dfdistri)) {
            scatterutil.star(coo\$x, coo\$y, dfdistri[, i], cstar = cstar,
                coul[i])
        }
    if (cellipse > 0)
        allmodel<-data.frame()
        for (i in 1:ncol(dfdistri)) {
              onecellipse<-data.frame()
              if (any(is.na(dfdistri[, i])))
                    return(invisible())
                if (sum(dfdistri[, i] * dfdistri[, i]) == 0)
                    return(invisible())

                coul = rep(1, length(coo\$x))
                cellipse = cellipse
                axesell = axesell

                util.ellipse <- function(mx, my, vx, cxy, vy, coeff) {
                lig <- 500
                epsi <- 1e-10
                x <- 0
                y <- 0
                if (vx < 0)
                    vx <- 0
                if (vy < 0)
                    vy <- 0
                if (vx == 0 && vy == 0)
                    return(NULL)
                delta <- (vx - vy) * (vx - vy) + 4 * cxy * cxy
                delta <- sqrt(delta)
                l1 <- (vx + vy + delta)/2
                l2 <- vx + vy - l1
                if (l1 < 0)
                    l1 <- 0
                if (l2 < 0)
                    l2 <- 0
                l1 <- sqrt(l1)
                l2 <- sqrt(l2)
                    test <- 0
                if (vx == 0) {
                    a0 <- 0
                    b0 <- 1
                    test <- 1
                    }
                if ((vy == 0) && (test == 0)) {
                    a0 <- 1
                    b0 <- 0
                    test <- 1
                }
        if (((abs(cxy)) < epsi) && (test == 0)) {
            a0 <- 1
            b0 <- 0
            test <- 1
        }
        if (test == 0) {
            a0 <- 1
            b0 <- (l1 * l1 - vx)/cxy
            norm <- sqrt(a0 * a0 + b0 * b0)
            a0 <- a0/norm
            b0 <- b0/norm
        }
        a1 <- 2 * pi/lig
        c11 <- coeff * a0 * l1
        c12 <- (-coeff) * b0 * l2
        c21 <- coeff * b0 * l1
        c22 <- coeff * a0 * l2
        angle <- 0
#        for (i in 1:lig) {
#            cosinus <- cos(angle)
#            sinus <- sin(angle)
#            x[i] <- mx + c11 * cosinus + c12 * sinus
#            y[i] <- my + c21 * cosinus + c22 * sinus
#            angle <- angle + a1
#        }
#        return(list(x = x, y = y, seg1 = c(mx + c11, my + c21,
#            mx - c11, my - c21), seg2 = c(mx + c12, my + c22,
#            mx - c12, my - c22)))
#         dx <- coeff * l1
#         dy <- coeff * l2
        return (c(mx,c11,c12,my,c21,c22))      #modified 2017.12.06,返回公式系数}
#         deta <- asin(b0)
#         return (c(mx,my,dx,dy,deta)) #modified 2017.12.21,返回中心点x,y坐标，长轴，短轴，旋转角度(正为逆时针，负为顺时针）}
    }
    dfdistri[, i] <- dfdistri[, i]/sum(dfdistri[, i])
    m1 <- sum(coo\$x * dfdistri[, i])
    m2 <- sum(coo\$y * dfdistri[, i])
    v1 <- sum((coo\$x - m1) * (coo\$x - m1) * dfdistri[, i])
    v2 <- sum((coo\$y - m2) * (coo\$y - m2) * dfdistri[, i])
    cxy <- sum((coo\$x - m1) * (coo\$y - m2) * dfdistri[, i])
    ell <- util.ellipse(m1, m2, v1, cxy, v2, cellipse)
    if (is.null(ell)){allmodel<-allmodel
    }else{
        #return(invisible())
#{    polygon(ell\$x, ell\$y, border = coul)
    #if (axesell)
        #segments(ell\$seg1[1], ell\$seg1[2], ell\$seg1[3], ell\$seg1[4],
            #lty = 2, col = coul)    #modified 2017.12.06}
    if (is.null(axesell)){
        allmodel<-rbind(allmodel,onemodel)
        }else{
# {       segments(ell\$seg2[1], ell\$seg2[2], ell\$seg2[3], ell\$seg2[4],
#            lty = 2, col = coul)
#            k<-length(ell\$x)
#            kfactor<-rep(i,k)
            kfactor <- i
            var = paste(ell,collapse=\",\")
            onemodel<-data.frame(kfactor,var)
            allmodel<-rbind(allmodel,onemodel)   #modified 2017.12.06}
            }
        }}
    if (clabel > 0)
        scatterutil.eti(coox, cooy, unique(as.factor(b)), clabel, coul = col)
    box()
    invisible(match.call())

    
    ##########
file1 <- paste(\"$opts{o}/\",paste(select,\"_circle.txt\",sep=\"\"),sep=\"\")
write.table(t(allmodel),file=file1,sep=\"\t\",quote = F,row.names=T,col.names=F)
#modified 2017.12.06

coo_zp<-data.frame(coo\$x,coo\$y)
dexy_zp<-data.frame(rownames(dfxy),dfxy)

#file2 <- paste(\"$opts{o}/\",paste(select,\"_point.txt\",sep=\"\"),sep=\"\")
#write.table(dexy_zp,file=file2,sep=\"\t\",quote = F,row.names=F,col.names=F)
#write.table(cbind(rownames(my_colum),my_colum),file=file2,sep=\"\t\",quote = F,row.names=F,col.names=c(\"Sample_ID\",\"x\",\"y\"))
#write.table(dexy_zp,file=file2,sep=\"\t\",quote = F,row.names=F,col.names=c(\"Sample_ID\",\"x\",\"y\"))

legend(\"topright\",xpd=TRUE,legend=groups,pch=pch,bty=\"n\")
if(select == \"BCA\"){
text(my_colum\$CS1+0.1,my_colum\$CS2+0.1,labels=names(data),font=2)
}else{
text(my_colum\$A1+0.1,my_colum\$A2+0.1,labels=names(data),font=2)
}
dev.off()

}

";





