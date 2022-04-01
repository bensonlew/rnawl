#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
my %opts;
my $VERSION = "v2.20171031"; 

GetOptions( \%opts,"i=s","o=s","c=s","i1=s","method=s","i2=s","pcnum=s");
my $usage = <<"USAGE";
       Program : $0   
       Discription:   
       Version : $VERSION
       Contact : hao.gao\@majorbio.com
       Usage :perl $0 [options]		
	   -i	*Input a table containing two column to do regression analysis	
           -o	Output dir
           -method data processing method;eg:center,scale,none
           -c	whether combine two table as input,default:F
           -i1	tax table
           -i2	function gene table
           -pcnum  The maximum number of calculations,default:3
USAGE

die $usage if(!($opts{i}||($opts{c}&&$opts{i1}&&$opts{i2})));
$opts{o}=defined $opts{o}?$opts{o}:"result";
$opts{c}=defined $opts{c}?$opts{c}:"F";
$opts{method}=defined $opts{method}?$opts{method}:"none";
$opts{pcnum}=defined $opts{pcnum}?$opts{pcnum}:"3";

open CMD,">RegressionAnalysis.cmd.r";
print CMD "

library(vegan)

mianfun<-function(diversity,pcLabel){
data<-data.frame(diversity)
#######中心化
if(\"$opts{method}\" == \"center\"){
data2<-scale(data,center=T,scale=F)
  }
####标准化
if(\"$opts{method}\" ==\"scale\"){
data2<-scale(data,center=T,scale=T)
}
####不处理
if(\"$opts{method}\" == \"none\"){
data2<-data
}
pcaone<-as.vector(data2\$diversity_tax)  # 因子变向量
pcaone<-as.numeric(pcaone)
gene<-as.vector(data2\$diversity_gene)  # 因子变向量
gene<-as.numeric(gene)
model<-lm(gene ~ pcaone)  # by houshuang, 横纵坐标对调
arry=summary(model)
data2<-data.frame(data2)
k<-arry\$coefficients[2]
b<-arry\$coefficients[1]
pcaonexmin=min(pcaone)
pcaonexmax=max(pcaone)
pcaoneymin=k*pcaonexmin+b
pcaoneymax=k*pcaonexmax+b
# by houshuang 20191009 >>>
adj_r_squ<-arry\$adj.r.squared
r_squ<-arry\$r.squared
fvalue<-arry\$fstatistic[1]
temp<-arry\$fstatistic
pvalue<-1-pf(temp[1],temp[2],temp[3])
# <<<
if(\"$opts{c}\" == \"F\"){
data2=data.frame(name=row.names(data2),func=data2\$diversity_gene,tax=data2\$diversity_tax)
}else{
data2=data.frame(name=data2\$sample,func=data2\$diversity_gene,tax=data2\$diversity_tax)
}
r_squ=data.frame(F=fvalue,p_value=pvalue,R.Square=r_squ,Adj.R.Square =adj_r_squ,xmin=pcaonexmin,ymin=pcaoneymin,xmax=pcaonexmax,ymax=pcaoneymax,K=k,B=b)
outfile1<-paste(\"$opts{o}.data\",pcLabel,\".xls\",sep='.')
outfile2<-paste(\"$opts{o}.message\",pcLabel,\".xls\",sep='.')
write.table(data2,outfile1,sep='\t',quote=F,row.names=F)
write.table(r_squ,outfile2,sep='\t',quote=F,row.names=F)
}

#######read diversity table
if(\"$opts{c}\" == \"F\"){
	diversity<-read.table(\"$opts{i}\",sep=\"\\t\",quote = \"\\t\",row.names=1,head=T,check.names = F)
	colnames(diversity)<- c(\"diversity_tax\",\"diversity_gene\")
    mianfun(diversity,\"PC\")

}else{
	tax<-read.table(\"$opts{i1}\",sep=\"\\t\",quote = \"\\t\",head=T,check.names = F)
	fun<-read.table(\"$opts{i2}\",sep=\"\\t\",quote = \"\\t\",head=T,check.names = F)
    tax<-tax[,1:2]
    colnames(tax)<-c(\"sample\",\"diversity_tax\")
    last_index = $opts{pcnum}+1
    if(length(fun[1,])<last_index){last_index=length(fun[1,])}
    for(i in seq(2,last_index,1)){
        fun_tmp<-fun[,c(1,i)]
        colnames(fun_tmp)<-c(\"sample\",\"diversity_gene\")
        diversity<-merge(tax,fun_tmp,by=\"sample\")
        mianfun(diversity,paste(\"PC\",i-1,sep=\"\"))
    system(\"cat $opts{o}.data.*.xls >  $opts{o}.data.xls\")
    system(\"cat $opts{o}.message.*.xls >  $opts{o}.message.xls\")
    #system(\"sed -i -e 's/tax/env/' -e 's/func/tax/' Regression.data.xls\") #guanqing 20180524 修改结果文件表头
    system(\"sed  -i 's/name\tfunc\ttax/name\ttax\tenv/' Regression.data.xls\")  #guanqing 20180525 修改结果文件表头

    }
}    

";

#`/mnt/ilustre/users/sanger-dev/app/program/R-3.3.1/bin/R  --restore --no-save < RegressionAnalysis.cmd.r`;
