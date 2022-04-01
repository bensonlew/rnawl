#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
my %opts;
my $VERSION = "V2.2017.3.22";
GetOptions( \%opts,"f=s","data=s","m=s","w=s","d=s","r=s","p=s","T=s","o=s","model=i","top=i","list=s");
my $usage = <<"USAGE";
       Program : $0   
       Discription:a multivariate statistical framework that finds associations between clinical metadata and microbial community abundance or function.   
       Version : $VERSION
       Contact : xuan.zou\@majorbio.com
       Usage :perl $0 -f env.txt -data otu_taxa_table.xls -m masslin.pcl -o Output
                        -f*             input clinical metadata,env.txt
                        -data*          input microbial community abundance or function,if model=0,input otu_taxa_table.xls.if model=1 input  tax or function table
                        -m*             out maaslin.pcl for the next analysis
                        -o*             the result of analysis
			-model		the model to analysis default:1,0:use all otus and taxes to analysis ;1:use one level to analysis
			-top		if model =1 ,you can select the top taxes to analysis ,default:100.
			-list		tax or function list,you can select the tax or function in which you are interested
                        -w              if model=0,you can select depth of taxonomy to be computed, negative=from right, 0=no change   default:-1
                        -d              dSignificanceLevel,anything equal to or lower than this is significant  default:0.05
                        -r              The minimum relative abundance allowed in the data. Values below this are removed and imputed as the median of the sample data  default:0.0001
                        -p              The minimum percentage of samples in which a feature must have the minimum relative abundance in order not to be removed     default:0.1
                        -G              This is the significance cuttoff of grubbs test that used to indicate an outlier or not. default:0.05
                        -T              testingCorrection,This indicates which multiple hypothesis testing method will be used, available are holm, hochberg, hommel, bonferroni, BH, BY.  default:BH
                        -l              transform_method [asinsqrt/none]  default:asinsqrt
                        -F              Metadata features that will be forced into the model seperated by commas. These features must be listed in the read.config file. Example '-F Metadata2,Metadata6,Metadata10'
                        -n              These data will not be imputed. Comma delimited data feature names. Example '-n Feature1,Feature4,Feature6'
                        
                        option with * must be seted.
                        Example:perl maaslin.pl -f env.txt -data otu_taxa_table.xls -m maaslin.pcl -o output

USAGE

die $usage if(!($opts{f}&&$opts{data}&&$opts{m}&&$opts{o}));
$opts{w}=defined $opts{w}?$opts{w}:-1;
$opts{d}=defined $opts{d}?$opts{d}:0.05;
$opts{r}=defined $opts{r}?$opts{r}:0.0001;
$opts{p}=defined $opts{p}?$opts{p}:0.1;
$opts{G}=defined $opts{G}?$opts{G}:0.05;
$opts{T}=defined $opts{T}?$opts{T}:"BH";
$opts{l}=defined $opts{l}?$opts{l}:"asinsqrt";
$opts{F}=defined $opts{F}?$opts{F}:"\"\"";
$opts{n}=defined $opts{n}?$opts{n}:"\"\"";
$opts{model}=defined $opts{model}?$opts{model}:1;
$opts{top}=defined $opts{top}?$opts{top}:100;
$opts{list}=defined $opts{list}?$opts{list}:"none";

#if ($opts{model}==0){
#    pass
#	#`python /mnt/ilustre/users/xuan.zou/test/qiimeToMaaslintest1.py $opts{f} <$opts{data}> $opts{m} -t $opts{w}`
#}
# add **fileEncoding='UTF-8'** in function read.table for unicode input
open RCMD, ">$opts{o}.cmd.r";
print RCMD "
library(Maaslin)
if ($opts{model}==1){
	top<-$opts{top}
	l<-\"$opts{list}\"
	env<-read.table(\"$opts{f}\",header=T,check.names=FALSE,sep=\"\\t\",fill=T,fileEncoding='UTF-8')
	rownames(env)<-as.character(env[,1])
        env<-env[,-1,drop=F]
	spe<-read.table(\"$opts{data}\",header=T,check.names=FALSE,sep=\"\\t\",fill=T,fileEncoding='UTF-8')
	rownames(spe)<-as.character(spe[,1])
	spe<-spe[,-1,drop=F]
	if(top >0){
		rsum <-sapply(1:nrow(spe),function(x) sum(spe[x,],na.rm=T))
		spe<-spe[order(rsum,decreasing=TRUE),]
		if($opts{top}<=nrow(spe)){
			spe<-spe[1:top,]
		}
	}
	if(l!=\"none\"){
		list_tax<-read.table(\"$opts{list}\",sep=\"\\t\",fileEncoding='UTF-8')
		spe<-spe[list_tax,]
	}
	abundance<-apply(spe,2,function(x) x/sum(x,na.rm=T))
	env<-t(env[colnames(abundance),,drop=F])
	pcl<-rbind(env,abundance)
	output_pcl<-rbind(colnames(pcl), pcl)
	rownames(output_pcl)[1]<-\"sample\"
	write.table(output_pcl,\"$opts{m}\",sep=\"\\t\",eol=\"\n\",quote=FALSE,col.names = F,row.names =T)
	writeLines(c(\"Matrix: Metadata\",paste(\"Read_PCL_Rows:-\",nrow(env)+1, sep = \"\", collapse = NULL),\"\",\"Matrix: Abundance\",paste(\"Read_PCL_Rows:\",nrow(env)+2,\"- \",sep = \"\", collapse = NULL)),\"maaslin_read.config\")
	}
Maaslin('$opts{m}','$opts{o}',strInputConfig=\"maaslin_read.config\",dSignificanceLevel=$opts{d},dMinAbd=$opts{r},dMinSamp=$opts{p},dPOutlier=$opts{G},strMultTestCorrection=\"$opts{T}\",strTransform=\"$opts{l}\",strForcedPredictors=$opts{F},strNoImpute=$opts{n})

#$`Rscript Maaslin.R  $opts{m} $opts{o} -i maaslin_read.config -d $opts{d} -r $opts{r} -p $opts{p} -G $opts{G} -T $opts{T} -l $opts{l} -F $opts{F} -n $opts{n}`;

######`R --restore --no-save < $opts{o}.cmd.r`;


##add line  20191010

mainfun<-function(diversity,pcLabel,spe){
data<-data.frame(diversity)

data2<-data

# diversity_tax : abund    diversity_gene: env
model<-lm(diversity_abund ~ diversity_env,data.frame(data2))
arry=summary(model)
data2<-data.frame(data2)
pcaone<-data2\$diversity_env
k<-arry\$coefficients[2]
b<-arry\$coefficients[1]
pcaonexmin=min(pcaone)
pcaonexmax=max(pcaone)
pcaoneymin=k*pcaonexmin+b
pcaoneymax=k*pcaonexmax+b
r_squ<-summary(model)\$adj.r.squared
#if(\"$opts{c}\" == \"F\"){
#data2=data.frame(name=row.names(data2),func=data2\$diversity_gene,tax=data2\$diversity_tax)
#}else{
#data2=data.frame(name=data2\$sample,func=data2\$diversity_gene,tax=data2\$diversity_tax)
#}
r_squ=data.frame(spe=spe, xmin=pcaonexmin,ymin=pcaoneymin,xmax=pcaonexmax,ymax=pcaoneymax,K=k,b=b)  #Adj.R.Square =r_squ,
return(r_squ)
#outfile1<-paste(\"$opts{o}.data\",pcLabel,\".xls\",sep='.')
#outfile2<-paste(\"$opts{o}.message\",pcLabel,\".xls\",sep='.')
#write.table(data2,outfile1,sep='\t',quote=F,row.names=F)
#write.table(r_squ,outfile2,sep='\t',quote=F,row.names=F)
}


#tax<-read.table(\"$opts{i1}\",sep=\"\\t\",quote = \"\\t\",head=T,check.names = F)
fun<-read.table(\"Maaslin/maaslin.tsv\",sep=\"\\t\",head=T,check.names = F,fileEncoding='UTF-8')
head_name<- names(fun)
tax<-fun[,c(1,2)]
colnames(tax)<-c(\"sample\",\"diversity_env\")
last_index = length(fun[1,])
#if(length(fun[1,])<last_index){last_index=length(fun[1,])}
all_result = 0
for(i in seq(3,last_index,1)){
	fun_tmp<-fun[,c(1,i)]
	spe <- head_name[i]
	colnames(fun_tmp)<-c(\"sample\",\"diversity_abund\")
	diversity<-merge(tax,fun_tmp,by=\"sample\")
	r_sum = mainfun(diversity,paste(\"spe\",i-1,sep=\"\"),spe)
	if(all_result == 0){
		all_result = r_sum
	}else{
		all_result = rbind(all_result, r_sum)
	}

}

outfile2<-paste(\"$opts{o}.message\",\"xls\",sep='.')
write.table(all_result,outfile2,sep='\t',quote=F,row.names=F)
#system(\"cat Maaslin/$opts{o}.message.*.xls >  $opts{o}.message.xls\")

#system(\"sed -i -e 's/tax/env/' -e 's/func/tax/' Regression.data.xls\") #guanqing 20180524 修改结果文件表头
#system(\"sed  -i 's/name\tfunc\ttax/name\ttax\tenv/' Regression.data.xls\")  #guanqing 20180525 修改结果文件表头

";
