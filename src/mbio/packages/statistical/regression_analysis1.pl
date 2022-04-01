#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
my %opts;
my $VERSION = "v2.20171031"; 

GetOptions( \%opts,"i=s","o=s","c=s","i1=s","method=s","i2=s");
my $usage = <<"USAGE";
       Program : $0   
       Discription:   
       Version : $VERSION
       Contact : hao.gao\@majorbio.com
       Usage :perl $0 [options]		
	   -i	*Input a table containing two column to do regression analysis	
           -o	Output dir
           -scale data processing method
           -c	whether combine two table as input,default:F
           -i1	tax table
           -i2	function gene table
USAGE

die $usage if(!($opts{i}||($opts{c}&&$opts{i1}&&$opts{i2})));
$opts{o}=defined $opts{o}?$opts{o}:"result";
$opts{c}=defined $opts{c}?$opts{c}:"F";
$opts{scale}=defined $opts{scale}?$opts{scale}:"F";

open CMD,">cmd.r";
print CMD "

library(vegan)

#######read diversity table
if(\"$opts{c}\" == \"F\"){
	diversity<-read.table(\"$opts{i}\",sep=\"\\t\",quote = \"\\t\",row.names=1,head=T,check.names = F)
	colnames(diversity)<- c(\"diversity_tax\",\"diversity_gene\")
}else{
	tax<-read.table(\"$opts{i1}\",sep=\"\\t\",quote = \"\\t\",head=T,check.names = F)
	fun<-read.table(\"$opts{i2}\",sep=\"\\t\",quote = \"\\t\",head=T,check.names = F)
	tax<-tax[,1:2]
	fun<-fun[,1:2]
	colnames(tax)<-c(\"sample\",\"diversity_tax\")
	colnames(fun)<-c(\"sample\",\"diversity_gene\")
	diversity<-merge(tax,fun,by=\"sample\")
}
data<-data.frame(diversity)

####标准化
if(\"$opts{scale}\" ==\"T\"){
data2<-scale(data,center=T,scale=T)
}
####不处理
if(\"$opts{scale}\" == \"F\"){
data2<-data
}
model<-lm(diversity_tax ~ diversity_gene,data.frame(data2))
arry=summary(model)
data2<-data.frame(data2)
pcaone<-data2\$diversity_gene
k<-arry\$coefficients[2]
b<-arry\$coefficients[1]
pcaonexmin=min(pcaone)
pcaonexmax=max(pcaone)
pcaoneymin=k*pcaonexmin+b
pcaoneymax=k*pcaonexmax+b
r_squ<-summary(model)\$adj.r.squared
if(\"$opts{c}\" == \"F\"){
data2=data.frame(Name=row.names(data2),x=data2\$diversity_gene,y=data2\$diversity_tax)
}else{
data2=data.frame(name=data2\$sample,func=data2\$diversity_gene,tax=data2\$diversity_tax)
}
r_squ=data.frame(Adj.R.Square =r_squ,xmin=pcaonexmin,ymin=pcaoneymin,xmax=pcaonexmax,ymax=pcaoneymax,K=k,b=b)
write.table(data2,\"$opts{o}_data.xls\",sep='\t',quote=F,row.names=F)
write.table(r_squ,\"$opts{o}_message.xls\",sep='\t',quote=F,row.names=F)
";

#`/mnt/ilustre/users/sanger-dev/app/program/R-3.3.1/bin/R  --restore --no-save < RegressionAnalysis.cmd.r`;
