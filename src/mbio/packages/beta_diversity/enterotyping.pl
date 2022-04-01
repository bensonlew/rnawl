#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
my %opts;
my $VERSION = "V1_20171027";

GetOptions( \%opts,"i=s","o=s","d=s");
my $usage = <<"USAGE";
       Program : $0   
       Discription:   
       Version : $VERSION
       Contact : hao.gao\@majorbio.com
       Last_motified : 2017.11.09 by zouxuan
       Usage :perl $0 [options]		
			-i	* input genus table file 
			-o	* output dir
			-d   dist_matrix

	   Example:$0 

USAGE
die $usage if(!($opts{i}&&$opts{o}));
$opts{d}=defined $opts{d}?$opts{d}:"none"; #add by zouxuan 20171120


### plot 
open RCMD, ">Enterotyping.r";
print RCMD "

library(\"cluster\")
library(\"clusterSim\")
library(\"ade4\")
data=read.table(\"$opts{i}\",sep=\"\\t\",header=T, row.names=1, check.names=F)
#data=data[-1,]
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
if(\"$opts{d}\" == \"none\"){
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
    nclusters[k]=0.00
  } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k]=index.G1(t(data),data.cluster_temp, d = data.dist, centrotypes = \"medoids\")
  }
}

k_clusters = c(1, 2, 3, 4, 5, 6, 7, 8, 9)
CH_index = nclusters
d <- cbind(k_clusters,CH_index)
file2 <- paste(\"$opts{o}/\",\"ch.txt\",sep=\"\")  #zx
write.table(d,file2,row.names=F,sep = \"\t\",quote = F) #zx
#print(head(nclusters))
#write.table(nclusters,file = \"ch.txt\", append = T,row.names=F,col.names=F,sep = \"\t\")
#m<-nclusters
#write.table(m,file = \".ch.txt\",append = T,quote = F,sep = \"\t\",eol = \"\t\",row.names=F,col.names=F) # zhouxuan

max=which.max(nclusters)
data.cluster=pam.clustering(data.dist, k=max)
Sample=names(data[0,])
Enterotype=data.cluster
all=cbind(Sample,Enterotype)

####�ó���?�������?###
file1 <- paste(\"$opts{o}/\",\"cluster.txt\",sep=\"\")
write.table(all,file1,sep=\"\t\",quote = F,row.names=F)

#####���ÿ���ȱ�######
for (i in 1:max){
  result=NULL
  s=NULL
 j<-all[,2] ==i
# if(length(all[j,])==2){
# all1<-data.frame(all)
# c<-all1[j,][1]
# c<-as.character(c)
# result<-data[,c]
# result\$sum<-result
# s=paste(round(result\$sum/sum(result\$sum)*100,2),\"%\",sep=\"\")
# result\$pre<-s
# #result=result[order(result\$sum,decreasing=T),]
# ss<- c(\"taxon_name\",names(result))
# file <- paste(\"$opts{o}/\",i,\".cluster.txt\",sep=\"\")
# #write.table(ss,file=file,sep=\"\t\",quote = F,eol=\"\t\",row.names=F,col.names=F)
# #write.table(\"\\r\",file=file,quote = F,sep=\"\t\",row.names=F,col.names=F,append=T)
# #write.table(result,file=file,sep=\"\t\",quote = F,row.names=T,append=T,col.names=F)
# }
# else{
 c<-all[j,,drop=F][,1]
 c<-as.character(c)
 result<-data[,c,drop=F]
 result\$sum<-rowSums(result)
 s=paste(round(result\$sum/sum(result\$sum)*100,2),\"%\",sep=\"\")
 result\$pre<-s
 result=result[order(result\$sum,decreasing=T),]
 ss<- c(\"taxon_name\",names(result))
 file <- paste(\"$opts{o}/\",i,\".cluster.txt\",sep=\"\")
 write.table(ss,file=file,sep=\"\t\",quote = F,eol=\"\t\",row.names=F,col.names=F)
 write.table(\"\\r\",file=file,quote = F,sep=\"\",row.names=F,col.names=F,append=T)
 write.table(result,file=file,sep=\"\t\",quote = F,row.names=T,append=T,col.names=F)
 }
 #c<-as.character(c) 
 #result<-data[,c]
 #result\$sum<-rowSums(result)
 #s=paste(round(result\$sum/sum(result\$sum)*100,2),\"%\",sep=\"\")
 #result\$pre<-s
 #result=result[order(result\$sum,decreasing=T),]
 #ss<- c(\"taxon_name\",names(result))
 #file <- paste(\"$opts{o}/\",i,\".cluster.txt\",sep=\"\")
 #write.table(ss,file=file,sep=\"\t\",quote = F,eol=\"\t\",row.names=F,col.names=F)
 #write.table(\"\\r\",file=file,quote = F,sep=\"\t\",row.names=F,col.names=F,append=T) 
 #write.table(result,file=file,sep=\"\t\",quote = F,row.names=T,append=T,col.names=F) 
#}



####��CH-indexָ��ͼ#####
#pdf=paste(\"CH-index\",\".pdf\",sep=\"\")
#pdf(pdf,width=8,height=5)
#par(mar=c(5,5,4,2))

#plot(nclusters, type=\"h\", xlab=\"k clusters\", ylab=\"CH index\",main=\"Optimal number of clusters\")
#dev.off()

";
#system ("/mnt/ilustre/users/sanger-dev/app/program/R-3.3.1/bin/Rscript Enterotyping.r");
#system ("/mnt/ilustre/users/sanger-dev/app/program/R-3.3.1/bin/R $opts{o}.Enterotyping.r");
#system ("Rscript $opts{o}.Enterotyping.r");


