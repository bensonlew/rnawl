#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my %opts;
my $VERSION = "v3.20180925";

GetOptions( \%opts,"mode=s","i=s","group=s","o=s","n=s","method=s","name=s","smooth=s","conf_level=s","print_auc=s","print_thres=s","labels=s","labelsize=s","ci=s","ci_se=s","ci_sp=s","h=f","w=s","type=s");
my $usage = <<"USAGE";
       Program : $0   
       Discription:pROC-package; Tools for visualizing, smoothing and comparing receiver operating characteristic (ROC curves). (Partial) area under the curve (AUC) can be compared with statistical tests based on U-statistics or bootstrap. Confidence intervals can be computed for (p)AUC or ROC curves.
       Version : $VERSION
       Contact : chengchen.ye\@majorbio.com
              Last_author : hao.gao
       Usage :perl $0 [options]         
                        -mode   *1 or 2 ; The Receiver:1)The relative abundance of the top n  Bacteria.
                                                       2)The relative abundance of special bacteria(one or more).    
                        -i      *Input genus table file(or other Taxonomic level).
                        -group  *Group name in mapping file. 1= case group, 0= control group.
                                the value of the response for controls and cases respectively. 
						(group.txt e.g.. samples        group1
                                                                 Sample1        1
                                                                 Sample2        0
                                                                 Sample3        1
                                                                 Sample4        0
                                                                  ...)

                        -o      *Output dir
                        -n      (*mode_1)Top n genus or other taxonomic level(relative abundance),Defult:50
                        -method (*mode_1)Follow method are available:sum, average, median, MI
                        -name   (*mode_2)the name of Bacteria 
                                         (name.txt e.g..name
                                                        Nitrospira
                                                        Acidibacter
                                                        Gemmatimonadaceae_uncultured
                                                        ...)
                        -conf_level     Confidence interval;(0,1) Default:0.97
			-print_auc	boolean. Should the numeric value of AUC be printed on the plot? Default:T
                        -print_thres	TRUE,FALSE or "best". Should a selected set of thresholds be displayed on the ROC curve?  Note that on a smoothed ROC curve, only "best” is supported.Default:"best"	
                        -ci	boolean. Should we plot the confidence intervals(best thresholds)? Defult:T 								
			-ci_se	boolean. This function adds "se" confidence intervals to a ROC curve plot.Defult:FALSE	
                        -ci_sp  boolean. This function adds "sp" confidence intervals to a ROC curve plot.Defult:FALSE 
                        -type	type of plot ci_se or ci_sp, “bars” or “shape”. Can be shortened to “b” or “s”. “shape” is only available for ci.se and ci_se and ci_sp, not for ci.thresholds.Defult:shape
                        -smooth	If smooth=TRUE,ci must be "FALSE". This function smoothes a ROC curve of numeric predictor. Default:FALSE.
                 	-w      default:7
                        -h      default:7														

           Example:perl roc_plus.pl -mode 1 -i otu_table.xls -group group.txt -o outdata -n 10 -method median  
                   perl roc_plus.pl -mode 2 -i otu_table.xls -group group.txt -o outdata -name name.txt 
USAGE

die $usage if(!($opts{mode}&&$opts{i}&&$opts{group}&&$opts{o}));

$opts{print_auc}=defined $opts{print_auc}?$opts{print_auc}:"TRUE";
$opts{print_thres}=defined $opts{print_thres}?$opts{print_thres}:"\"best\"";
$opts{smooth}=defined $opts{smooth}?$opts{smooth}:"FALSE";
$opts{ci}=defined $opts{ci}?$opts{ci}:"TRUE";
$opts{ci_se}=defined $opts{ci_se}?$opts{ci_se}:"FALSE";
$opts{ci_sp}=defined $opts{ci_sp}?$opts{ci_sp}:"FALSE";
$opts{type}=defined $opts{type}?$opts{type}:"shape";
$opts{n}=defined $opts{n}?$opts{n}:"20";
$opts{method}=defined $opts{method}?$opts{method}:"sum";
$opts{name}=defined $opts{name}?$opts{name}:"NULL";
$opts{ncuts}=defined $opts{ncuts}?$opts{ncuts}:20;
$opts{labels}=defined $opts{labels}?$opts{labels}:"F";
$opts{labelsize}=defined $opts{labelsize}?$opts{labelsize}:"3";
$opts{w}=defined $opts{w}?$opts{w}:7;
$opts{h}=defined $opts{h}?$opts{h}:7;
$opts{conf_level}=defined $opts{conf_level}?$opts{conf_level}:"0.97";

if(! -e $opts{o}){
                `mkdir $opts{o}`;
}


open CMD,">$opts{o}/roc.cmd.r";

print CMD "
library(pROC)
##############################################################################################
##############################################################################################
if($opts{mode}==1){
if(\"$opts{method}\"==\"MI\" ){
group <- read.delim(\"$opts{group}\",header=T,check.names=F,comment.char=\"\",colClasses=c('character'))
otu_table <-read.delim(\"$opts{i}\",sep=\"\\t\",head=T,check.names = F,row.names=1)
rocdata<-as.data.frame(group)
rocdata[,2]<- as.numeric(rocdata[,2])
a<-length(otu_table[,1])
b<-$opts{n}
if(b==0){
b<-length(otu_table[,1])
}
c<-min(a,b)
}
else{
group <- read.delim(\"$opts{group}\",header=T,check.names=F,comment.char=\"\",colClasses=c('character'))
otu_table <-read.delim(\"$opts{i}\",sep=\"\\t\",head=T,check.names = F,row.names=1)
###choose samples to be analyzed throuth group.txt
library(dplyr)
otu_table2<-as.data.frame(otu_table)
otu_table2[is.na(otu_table2)]<-0
###select top n 
rowsum <-sapply(1:nrow(otu_table2),function(x) sum(otu_table2[x,]))
otu_table3<-otu_table2[order(rowsum,decreasing=TRUE),]
a<-length(otu_table[,1])
b<-$opts{n}
if(b==0){
b<-length(otu_table[,1])
}
c<-min(a,b)
top_n<-otu_table3[1:c,]
###############################method=sum,average,median
top_n<-t(top_n)
############sum 
if(\"$opts{method}\"==\"sum\" ){
datasum<-sapply(1:nrow(top_n),function(x) sum(top_n[x,]))
}
############average
if(\"$opts{method}\"==\"average\" ){
#datasum<-sapply(1:nrow(top_n),function(x) sum(top_n[x,])/$opts{n})
datasum<-sapply(1:nrow(top_n),function(x) sum(top_n[x,])/c)
}
############median 
if(\"$opts{method}\"==\"median\" ){
datasum<-sapply(1:nrow(top_n),function(x) median(top_n[x,]))
}
############end
datasum<-as.data.frame(datasum)
datasum<-cbind(as.data.frame(row.names(top_n)),datasum)
colnames(datasum)<-c(\"samples\",\"sum\")
names(group)<-c(\"samples\",\"group1\")
############################end

rocdata<-merge(datasum,group,by=\"samples\",all=T)
rocdata<-as.data.frame(rocdata)

#############################plot
rocdata[,2]<- as.numeric(rocdata[,2])
}

otu<-plot.roc(rocdata[,3], rocdata[,2],direction=c(\"auto\"),ci=T,of=\"thresholds\",thresholds = \"best\",specificities = seq(0, 100, 5),percent =T)
roc_normal_otu2 <- function(){
    otu2 <<- plot.roc(rocdata[,3], rocdata[,2],ci=T,percent =T,smooth=T)
}
tryCatch(roc_normal_otu2(), error=function(e){print(\"error info\"); print(e)})

####画曲线点并输出结果
#roctable<-cbind(otu\$specificities,otu\$sensitivities)
roc_normal_run <- function(){
    roctable <<- cbind(otu\$specificities,otu\$sensitivities)
    colnames(roctable) <<- c(\"specificity\",\"sensitivity\")
    table <- paste(\"$opts{o}/\",\"roc_curve.xls\",sep=\"\")
    write.table(roctable, table, row.names=F,sep=\"\\t\",quote =F)
    }
tryCatch(roc_normal_run(), error=function(e){print(\"error info\"); print(e)})


###平滑处理曲线点
#roctable2<-cbind(otu2\$specificities,otu2\$sensitivities)
roc_normal_run2 <- function(){
    roctable2 <<- cbind(otu2\$specificities,otu2\$sensitivities)
    colnames(roctable2)<-c(\"specificity\",\"sensitivity\")
    table1 <- paste(\"$opts{o}/\",\"roc_curve_smooth.xls\",sep=\"\")
    write.table(roctable2, table1, row.names=F,sep=\"\\t\",quote =F)
    }
tryCatch(roc_normal_run2(), error=function(e){print(\"error info\"); print(e)})

#####输出置信区间
#rocarea_value<-ci.se(otu,specificities = seq(0, 100, 2),boot.n = 1000,conf.level=$opts{conf_level}, stratified=FALSE)
rocarea_normal_run <- function(){
    rocarea_value <- ci.se(otu,specificities = seq(0, 100, 2),boot.n = 1000,conf.level=$opts{conf_level}, stratified=FALSE)
    c1<-(1-$opts{conf_level})/2*100
    c2<- \"50%\"
    c3<-(1-(1-$opts{conf_level})/2)*100
    cc1<-paste(\'sensitivity\',\'.\',c1,\'%\',sep = \"\")
    cc2<-paste(\'sensitivity\',\'.\',c2,sep = \"\")
    cc3<-paste(\'sensitivity\',\'.\',c3,\'%\',sep = \"\")
    curve <- paste(\"$opts{o}/\",\"roc_interval.xls\",sep=\"\")
    rocarea_value<-as.data.frame(rocarea_value)
    rocarea_value2<-data.frame(specificity=row.names(rocarea_value),aa=rocarea_value[,1],bb=rocarea_value[,2],cc=rocarea_value[,3])
    colnames(rocarea_value2)=c(\"specificity\",cc1,cc2,cc3)
    write.table(rocarea_value2, curve, row.names=FALSE,sep=\"\t\",quote =F)
    }
tryCatch(rocarea_normal_run(), error=function(e){print(\"error info\"); print(e)})

#####输出最佳点
best_normal_run <- function(){
    best_loc<-ci(otu,of=\"thresholds\", thresholds=\"best\",conf.level=$opts{conf_level})
    sp1<-(1-$opts{conf_level})/2*100
    sp2<-\"50%\"
    sp3<-(1-(1-$opts{conf_level})/2)*100
    se1<-(1-$opts{conf_level})/2*100
    se2<-\"50%\"
    se3<-(1-(1-$opts{conf_level})/2)*100
    cc1<-paste(\'specificity\',\'.\',sp1,\'%\',sep = \"\")
    cc2<-paste(\'specificity\',\'.\',sp2,sep = \"\")
    cc3<-paste(\'specificity\',\'.\',sp3,\'%\',sep = \"\")
    cc4<-paste(\'sensitivity\',\'.\',se1,\'%\',sep = \"\")
    cc5<-paste(\'sensitivity\',\'.\',se2,sep = \"\")
    cc6<-paste(\'sensitivity\',\'.\',se3,\'%\',sep = \"\")
    best_loc<-as.data.frame(best_loc)
    best_loc1 <- data.frame(cut_off=row.names(best_loc),cc1=best_loc[1,1],cc2=best_loc[1,2],cc1=best_loc[1,3],cc1=best_loc[1,4],cc1=best_loc[1,5],cc1=best_loc[1,6])
    colnames(best_loc1)=c(\"specificity\",cc1,cc2,cc3,cc4,cc5,cc6)
    bestloc <- paste(\"$opts{o}/\",\"best_loc.xls\",sep=\"\")
    write.table(best_loc1,bestloc,row.names=FALSE,sep=\"\\t\",quote =F)
}
tryCatch(best_normal_run(), error=function(e){print(\"error info\"); print(e)})

#####输出AUC值
auc_normal_run <- function(){
    ci_auc<-ci.auc(otu,conf.level=$opts{conf_level})
    ci_auc<-as.data.frame(ci_auc)
    auc_data<-data.frame(conf_level=$opts{conf_level},min_auc=ci_auc[1,1],mean_auc=ci_auc[2,1],max_auc=ci_auc[3,1])
    auc<- paste(\"$opts{o}/\",\"roc_auc.xls\",sep=\"\")
    write.table(auc_data, auc, row.names=FALSE,sep=\"\\t\",quote =F)
}
tryCatch(auc_normal_run(), error=function(e){print(\"error info\"); print(e)})

#####输出piAUC值
pi_auc_run <- function(){
    ci_auc2<-ci.auc(otu2,conf.level=$opts{conf_level})
    ci_auc2<-as.data.frame(ci_auc2)
    auc_data2<-data.frame(conf_level=$opts{conf_level},min_auc=ci_auc2[1,1],mean_auc=ci_auc2[2,1],max_auc=ci_auc2[3,1])
    auc2<- paste(\"$opts{o}/\",\"roc_auc_smooth.xls\",sep=\"\")
    write.table(auc_data2, auc2, row.names=FALSE,sep=\"\\t\",quote =F)
}
tryCatch(pi_auc_run(), error=function(e){print(\"error info\"); print(e)})
}

";


