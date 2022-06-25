#!/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/perl
use strict;
use warnings;
use Getopt::Long;
my %opts;
my $VERSION = "V2.20160708";
GetOptions( \%opts,"i=s","m=s","o=s","g=s","ntree=i","top=i","type=s");
my $usage = <<"USAGE";
       Program : $0   
       Discription: Program used to caculate randomforest,with mds plot and importance variables given .  
       Version : $VERSION
       Contact : jie.yao\@majorbio.com  # last modify by shaohua.yuan in 20171103
       Usage :perl $0 [options]         
                        -i      * input otu table file 
                        -o      * output dir
                        -m      input mapping file if you want set points\'s color and pch by groups. If omitted, randomForest will run in unsupervised mode.
				Default:none
                        -g      group name in mapping file .Default:none
                        -ntree	Number of trees to grow. This should not be set to too 	small a number, to ensure that every input row gets predicted at least a few times.Default:500	
        	            -top    How many variables to show? 
                        -type 	either 1,2 or 3, specifying the type of importance measure (1=mean decrease in accuracy, 2=mean decrease in node impurity).
                        -scal   importance should be scaled. Default: T. # add by shaohua.yuan
			
       Example:$0 -i otu_table.xls -o randomForest  -m group -g group  

USAGE

die $usage if(!($opts{i}&&$opts{o}));
die $usage if($opts{m}&& !$opts{g});
die $usage if(!$opts{m}&& $opts{g});

$opts{m}=defined $opts{m}?$opts{m}:"none";
$opts{g}=defined $opts{g}?$opts{g}:"none";
$opts{ntree}=defined $opts{ntree}?$opts{ntree}:"500";
$opts{type}=defined $opts{type}?$opts{type}:"1";
$opts{top}=defined $opts{top}?$opts{top}:"50";
$opts{scal}=defined $opts{scal}?$opts{scal}:"T";

if(! -e $opts{o}){
                `mkdir $opts{o}`;
}


open CMD,">$opts{o}/cmd.r";
print CMD "
library(sp,warn.conflicts = F)
library(randomForest,warn.conflicts = F)
library(maptools,warn.conflicts = F)
library(stringr)
basename=\"randomforest\"


# if read otu data
otu <-read.table(\"$opts{i}\",sep=\"\\t\",head=T,check.names = F)
b1<-ncol(otu)
a1<-nrow(otu)
otu1<-matrix(nrow=a1,ncol=1)
#for(i in c(2:b1)){
#for(j in c(1:a1)){
#otu[j,i]<-otu1[j,i]
#}
#}


#otu3<-c()

zp<-length(str_split(otu[1,1],';')[[1]])

for(i in c(1:nrow(otu))){
otu1[i,1]<-str_split(otu[i,1],';')[[1]][zp]
#otu1[i,1]<-sub(':',otu1[i,1])
#otu1[i,1]<-str_split(otu[i,1],';')[[1]][2]
}

for(i in c(1:nrow(otu))){
#otu1[i,1]<-str_split(otu[i,1],':')[[1]][1]
otu1[i,1]<-sub(':','_',otu1[i,1])
#otu1[i,1]<-str_split(otu[i,1],';')[[1]][2]
}


otu\$otu_id<-otu1
otu<-otu[,-1]

abcd<-colnames(otu)
otu5<-data.frame(otu\$otu_id,otu[,-ncol(otu)])

names(otu5)<-c(abcd[length(abcd)],abcd[-length(abcd)])
#abcd<-colnames(otu)
otu<-otu5

rownames(otu) <-as.factor(otu[,1])
otu <-otu[,-1]
rownames(otu) <-sapply(rownames(otu),function(x) gsub(\"_*{.+}\",\" \",x,perl = TRUE))
rownames(otu) <-sapply(rownames(otu),function(x) gsub(\"-\",\"_\",x,perl = TRUE))
rownames(otu) <-sapply(rownames(otu),function(x) gsub(\"\\\\[\",\"\",x,perl = TRUE))
rownames(otu) <-sapply(rownames(otu),function(x) gsub(\"\\\\]\",\"\",x,perl = TRUE))
rownames(otu) <-sapply(rownames(otu),function(x) gsub(\"\\\\(\",\"\",x,perl = TRUE))
rownames(otu) <-sapply(rownames(otu),function(x) gsub(\"\\\\)\",\"\",x,perl = TRUE))
#rownames(otu) <-sapply(rownames(otu),function(x) gsub(\"^[0-9]\",\"X\\\\1\",x,perl = TRUE))
rownames(otu) <-sapply(rownames(otu),function(x) gsub(\"\/\",\"\",x,perl = TRUE))
otu <-as.data.frame(t(otu),stringsAsFactors=T)

map=\"$opts{m}\"
if(map !=\"none\"){
                sd <-read.table(\"$opts{m}\",head=T,sep=\"\\t\",comment.char = \"\",check.names = FALSE)        
                rownames(sd) <- as.character(sd[,1])
                sd[,1] <-as.character(sd[,1])
                sd\$group <-as.factor(sd[,2] )
                legend <- as.matrix(unique(sd\$group)) 
}

set.seed(1)
if(map != \"none\"){
otu.rf <- randomForest(sd\$group ~ .,otu,importance=T,proximity=T,ntree=$opts{ntree})



class_count <-as.matrix(table(sd\$group))
class <-data.frame(count=class_count)

##randomforest votes probably
votes_probably<- paste(\"$opts{o}/\",basename,\"_votes_probably.xls\",sep=\"\")
write.table(otu.rf\$votes,votes_probably,sep=\"\\t\",quote=F)

##randomforest predicted answer
predicted_answer <- paste(\"$opts{o}/\",basename,\"_predicted_answer.xls\",sep=\"\")
write.table(otu.rf\$predicted,predicted_answer,sep=\"\\t\",quote=F)

##randomforest classification table
rf_table <- paste(\"$opts{o}/\",basename,\"_confusion_table.xls\",sep=\"\")
write.table(otu.rf\$confusion,rf_table,sep=\"\\t\",quote=F)
mds <- cmdscale(1-otu.rf\$proximity)  
}else{
otu.rf <- randomForest(otu,importance=T,proximity=T,ntree=$opts{ntree})
mds <- cmdscale(1-otu.rf\$proximity)
}


##mds points
mds_points <- paste(\"$opts{o}/\",basename,\"_mds_sites.xls\",sep=\"\")
write.table(mds,mds_points,sep=\"\\t\",quote=F)

##proximity table
proximity <- paste(\"$opts{o}/\",basename,\"_proximity_table.xls\",sep=\"\")
write.table(otu.rf\$proximity,proximity,sep=\"\\t\",quote=F)

## importance table
vimp_table <- paste(\"$opts{o}/\",basename,\"_vimp_table.xls\",sep=\"\")
#write.table(otu.rf\$importancei,vimp_table,sep=\"\\t\",quote=F)
#imp<-importance(otu.rf)
#imp[,\"MeanDecreaseAccuracy\"]<-abs(imp[,\"MeanDecreaseAccuracy\"]) modify by shaohua.yuan
imp<-importance(otu.rf,scale = $opts{scal})
if($opts{type} == 1){
    imp <- imp[order(imp[,\"MeanDecreaseAccuracy\"],decreasing=T),]
    }else{
    imp <- imp[order(imp[,\"MeanDecreaseGini\"],decreasing=T),]
    }
write.table(imp,vimp_table,sep=\"\\t\",quote=F)


## top importance species table
topx_vimp <- paste(\"$opts{o}/\",basename,\"_topx_vimp.xls\",sep=\"\")
#imp <- importance(otu.rf)
imp <- importance(otu.rf,scale = $opts{scal})
if($opts{type} == 1){
	top <- imp[order(imp[,\"MeanDecreaseAccuracy\"],decreasing=T),][1:min($opts{top},length(imp[,1])),]
	write.table(t(otu)[rownames(top),],topx_vimp,sep=\"\\t\",quote=F)
	
}else if ($opts{type} == 2){
	top <- imp[order(imp[,\"MeanDecreaseGini\"],decreasing=T),][1:min($opts{top},length(imp[,1])),]
        write.table(t(otu)[rownames(top),],topx_vimp,sep=\"\\t\",quote=F)
}

";

`R --restore --no-save < $opts{o}/cmd.r`;
