#!/usr/bin/perl -W

use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions (\%opts,"spe=s","env=s","g=s","o=s","filter=s");
my $usage = <<"USAGE";
	Version :2017.01
	Contact : juan.zhu\@majorbio.com
	Usage	:perl $0 [options]
		-spe*	species table
		-env* 	environment table
		-g*	lable design file,only 2~4 group
			env1        soil
                        env2        soil
                        env3        soil
                        ENV1        air
                        ENV2        air
                        ENV3        air
		-o*	the name of pdf
		-filter filter env by group ,T or F ,default:F
	Example:perl VPA.pl -spe genus.xls -env env.xls -g group.xls -o test_genus

USAGE
die $usage if ( !($opts{spe}&&$opts{env}&&$opts{g}&&$opts{o}));
$opts{filter}=defined $opts{filter}?$opts{filter}:"F";


open CMD,">VPA.cmd.r";
print CMD "
spe <-read.table(file=\"$opts{spe}\",header=T,row.names=1,check.names=FALSE,sep=\"\\t\",quote=\"\",comment.char=\'\')
env <-read.table(file=\"$opts{env}\",header=T,row.names=1,check.names=FALSE,sep=\"\\t\",comment.char=\'\')
env <- env[which(rownames(env) %in% colnames(spe)),]
#####给环境因子分组，且只能是2~4组#####
group <- read.table(file=\"$opts{g}\",sep=\"\\t\")
glst <- lapply(1:length(unique(group[,2])),function(x)group[which(group[,2] %in% unique(group[,2])[x]),1])
names(glst) <-as.character(unique(group[,2]))
n<-length(names(glst))

#####根据env挑选样品#####
spe_last <- spe[,which (colnames(spe) %in% rownames(env))]

#####筛选环境因子#####
library(vegan)
filter<-\"$opts{filter}\"
if (filter == \"T\"){
for (i in 1:length(names(glst))){
	tmp <-env[,which (colnames(env) %in% glst[[i]])]
	stepforward<-ordistep(rda(t(spe_last)~ 1,data=tmp),scope=formula(rda(t(spe_last) ~ .,data=tmp)),direction=\"forward\",pstep=1000)
env.vif<-vif.cca(stepforward) #stepforward后剩下的环境因子
glst[[i]]<-row.names(as.data.frame(env.vif))
	}
}

#####计算每个环境因子的r.squared adj.r.squared#####
R2adj <- c()
for (i in  1:length(colnames(env))){
	spe.rda <-rda(t(spe_last)~env[,i])
	R2adj <- c(R2adj,RsquareAdj(spe.rda)\$adj.r.squared)
}

env_Radj <-t(rbind(colnames(env),R2adj))
env_Radj<-as.data.frame(env_Radj)
env_Radj <- merge(env_Radj,group,by.x=\'V1\',by.y=\'V1\')
colnames(env_Radj) <- c(\'Env\',\'Radj\',\'Group\')
write.table(env_Radj,\"env.R2adj.xls\",sep=\"\t\",col.names=T,quote=F,row.names = F)

#####增加图例#####
#env_Radj<- env_Radj[-1,]
#for (i in 1:length(names(glst))){
#       env_Radj[i]<- env_Radj[,which(colnames(env_Radj) %in% glst[[i]])]
#}

#####按组数作图#####
if(n==2){
#showvarparts(2)
mod <-varpart(t(spe_last),env[,which (colnames(env) %in% glst[[1]])],env[,which (colnames(env) %in% glst[[2]])])
mycolor=c(\"red\",\"blue\")
}
if(n==3){
#showvarparts(3)
mod <-varpart(t(spe_last),env[,which (colnames(env) %in% glst[[1]])],env[,which (colnames(env) %in% glst[[2]])],env[,which (colnames(env) %in% glst[[3]])])
mycolor=c(\"red\",\"blue\",\"green\")
}
if(n==4){
#showvarparts(4)
mod <-varpart(t(spe_last),env[,which (colnames(env) %in% glst[[1]])],env[,which (colnames(env) %in% glst[[2]])],env[,which (colnames(env) %in% glst[[3]])],env[,which (colnames(env) %in% glst[[4]])])
mycolor=c(\"red\",\"blue\",\"green\",\"yellow\")
}
pdf(\"$opts{o}.pdf\")
plot(mod, cutoff = -Inf,Xnames=names(glst),bg=mycolor)
dev.off()

#svg(\"$opts{o}.svg\")
#plot(mod, cutoff = -Inf,Xnames=names(glst),bg=mycolor)
#dev.off()


#write.table(mod\$part\$indfract\$Adj.R.square,\'env.plot.xls\',sep=\'\t\')  #guanqing add
env_plot = data.frame(lab=rownames(mod\$part\$indfract),Adj.R.squared=mod\$part\$indfract\$Adj.R.square)
write.table(env_plot,\'env.plot.xls\',sep=\'\t\',row.names =F, quote=F)

new_lab = \"group_name\"
for (each_lab in names(glst)){
    new_lab = paste(new_lab, each_lab, sep=\';\')
}

print(new_lab)

cmd = paste(\"sed -i \'1i \" , new_lab, \" \' env.plot.xls\", sep=\' \' )
system(cmd)

";
close CMD;
#`R --restore --no-save <VPA.cmd.r`;
