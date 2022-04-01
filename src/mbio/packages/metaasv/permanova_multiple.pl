#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
my $VERSION = "2017-2-20";
GetOptions (\%opts,"i=s","env=s","pe=i","o=s","rg=s","rtop=i","method=s","d_method=s","binary=s","multi=s");
my $usage = <<"USAGE";
Program : $0
Version : $VERSION
Contact : xuan.zou\@majorbio.com
Discription: Adnois / PERMANOVA analysis
Usage:perl $0 [options]
         data options:
         	-i *          input data matrix
         	-env *        env.txt/map.txt
         	-o *          output dir name
         	-rg	      row group,given a file with two column[tax1     group1	tax2     group2] [tax3	group3]
			-pe           permutation times,default : 999
			-rtop         choose the max top rows ,default:0
			-method	      use the method to get the q value ,default:BH
			-d_method     use the method to calculate distance,default:bray
        	-binary       if or not binary default :false
        	-multi		  output file name

USAGE
die $usage if ( !(defined $opts{i}&&$opts{env}&&$opts{o}&&$opts{multi}));
$opts{pe}=defined$opts{pe}?$opts{pe}:999;
$opts{rtop}=defined$opts{rtop}?$opts{rtop}:0;
$opts{rg}=defined$opts{rg}?$opts{rg}:"none";
$opts{method}=defined$opts{method}?$opts{method}:"BH";
$opts{d_method}=defined$opts{d_method}?$opts{d_method}:"bray";
$opts{binary}=defined$opts{binary}?$opts{binary}:"false";
my @lock=("GROUP","HEIGHT");
#my $env;
my $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   chop($script_dir);
   $script_dir = "./" unless ($script_dir);
open RCMD, ">cmd.r";
print RCMD "
library(vegan)
pe<-$opts{pe};rtop<-$opts{rtop}
spe<-read.table(file=\"$opts{i}\",header=T,check.names=FALSE,sep=\"\\t\",fill=T,comment.char=\"\")
rownames(spe)<-as.character(spe[,1])
spe<-spe[,-1,drop=F]
rgs <-\"$opts{rg}\"
if(rgs!=\"none\"){
	group <- read.table(\"$opts{rg}\")
      	glst <- lapply(1:length(unique(group[,2])),function(x)group[which(group[,2] \%in\% unique(group[,2])[x]),1])
      	names(glst) <-unique(group[,2])
      	tab <-t(sapply(1:length(glst),function(x) apply(spe[as.character(as.vector(glst[[x]])),,drop=F],2,mean)))
      	spe <-tab      
      	rownames(spe) <-unique(group[,2])
}else{
        if(rtop >0){
                rsum <-sapply(1:nrow(spe),function(x) sum(spe[x,]))
                spe<-spe[order(rsum,decreasing=TRUE),]
                if($opts{rtop}<=nrow(spe)){
                spe<-spe[1:rtop,]
        	}
	}
}
###单因素分析
env<-read.table(file=\"$opts{env}\",header=T,check.names=FALSE,sep=\"\\t\",fill=T,comment.char=\"\")
rownames(env)<-as.character(env[,1])
env<-env[,-1,drop=F]
spe<-spe[,rownames(env),drop=F]
spetmp<-t(spe)
p_value<-c()
R2<-c()
Df<-c()
SumsOfSqs<-c()
MeanSqs<-c()
F_Model<-c()
for (j in 1:ncol(env)){
        f<- adonis(spetmp ~ env[,j], data=env, permutations=pe,method=\"$opts{d_method}\",binary=\"$opts{binary}\")
        p_value<-c(p_value,f\$aov.tab\$Pr[1])
        R2<-c(R2,f\$aov.tab\$R2[1])
	Df<-c(Df,f\$aov.tab\$Df[1])
	SumsOfSqs<-c(SumsOfSqs,f\$aov.tab\$SumsOfSqs[1])
	MeanSqs<-c(MeanSqs,f\$aov.tab\$MeanSqs[1])
	F_Model<-c(F_Model,f\$aov.tab\$F.Model[1])
}
Characteristics<-colnames(env)
p_adjust=p.adjust(p_value,method=\"$opts{method}\",n=length(p_value))
all<-cbind(Characteristics,SumsOfSqs,MeanSqs,F_Model,R2,p_value,p_adjust)
PERMANOVA<-all[order(all[,6]),,drop=F]
write.table(PERMANOVA,\"$opts{o}\",sep=\"\t\",eol=\"\n\",quote=FALSE,col.names = T,row.names =F)

###多因素分析
env<-read.table(file=\"$opts{env}\",header=T,check.names=FALSE,sep=\"\\t\",fill=T,comment.char=\"\")
rownames(env)<-as.character(env[,1])
env<-env[,-1,drop=F]
spe<-spe[,rownames(env),drop=F]
spetmp<-t(spe)
x <- paste(colnames(env), collapse=\" + \")
method <- \"$opts{d_method}\"
binary <- \"$opts{binary}\"
fss <- paste(\"adonis(spetmp ~\", x, \",data=env, permutations=\", pe, \",method=\'\",method,\"\',binary=\'\", binary, \"\')\", sep=\"\")
f <- eval(parse(text=fss))
p_value<-c()
R2<-c()
Df<-c()
SumsOfSqs<-c()
MeanSqs<-c()
F_Model<-c()
p_value<-c(p_value,f\$aov.tab\$Pr)
R2<-c(R2,f\$aov.tab\$R2)
Df<-c(Df,f\$aov.tab\$Df)
SumsOfSqs<-c(SumsOfSqs,f\$aov.tab\$SumsOfSqs)
MeanSqs<-c(MeanSqs,f\$aov.tab\$MeanSqs)
F_Model<-c(F_Model,f\$aov.tab\$F.Model)
#Characteristics<-colnames(env)
#vector_name <- c(\"Residuals\", \"Total\")
#Characteristics <- c(Characteristics, vector_name)
p_adjust=p.adjust(p_value,method=\"$opts{method}\",n=length(p_value))
#all<-cbind(Characteristics,Df,SumsOfSqs,MeanSqs,F_Model,R2,p_value,p_adjust)
all<-cbind(Df,SumsOfSqs,MeanSqs,F_Model,R2,p_value,p_adjust)
MULTIPLE<-all[order(all[,6]),,drop=F]
write.table(MULTIPLE,\"$opts{multi}\",sep=\"\t\",eol=\"\n\",quote=FALSE,col.names = T,row.names =F)
write.table(f\$aov.tab,\"test_result.xls\",sep=\"\t\",eol=\"\n\")
";
