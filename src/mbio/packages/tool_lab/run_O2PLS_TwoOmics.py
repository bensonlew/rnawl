# -*- coding: utf-8 -*-

import re
import sys
import argparse
import os
from mako.template import Template


parser = argparse.ArgumentParser(description="run O2PLS analysis for two omics correlations")
parser.add_argument("-omics_x", type=str, required=True, help="the first omics data matrics")
parser.add_argument("-omics_y", type=str, required=True, help="the second omics data matrics")
parser.add_argument("-group",type =str, required=True, help = "group.txt")
parser.add_argument("-log_x",type =str, default="10", help = "log10,log1 no log")
parser.add_argument("-log_y",type =str, default="10", help = "log10,log1 no log") 
parser.add_argument("-oxoy", type=str, required=True, help="the two omics name such as 'Microbe;Metabolome'")
parser.add_argument("-scale", type=str, default="Par", help="the scale method to deal data 'UV',Ctr,Par,None")
parser.add_argument("-out", type=str, default="O2PLS", help="the output prefix")
parser.add_argument("-path", type=str, required=True, help="r file")
args = parser.parse_args()

Rcmd=r"""
library(OmicsPLS)
library(ggplot2)
library(dplyr)
library(snow)
source("${path}")

ScaleUnit<-function(x){
    x /sqrt(sum(x^2));
}


O2O <- "${oxoy}"
scalenm <- "${scale}"
O2O_sp <- strsplit(O2O,";")[[1]]
OX_dat <- read.table("${omics_x}",header = T,sep="\t",row.names=1)
OY_dat <- read.table("${omics_y}",header =T, sep="\t",row.names = 1,quote='\"')
group <- read.table("${group}",sep="\t")
names(group) <-c("Sample","Group")

if(as.numeric(${log_x}) >1){OX_dat = log(OX_dat+1,base=as.numeric(${log_x}))
}else{OX_dat = OX_dat}
if(as.numeric(${log_y}) > 1){
OY_dat = log(OY_dat+1,base=as.numeric(${log_y}))
}else{OY_dat = OY_dat}

centered_OX_data = t(scaleNorm(t(OX_dat),scalenm))
centered_OY_data = t(scaleNorm(t(OY_dat),scalenm))
#centered_OX_data = t(scale(t(OX_dat), scale=T))
#centered_OY_data = t(scale(t(OY_dat), scale=T))
ccc_OX_data <- centered_OX_data[complete.cases(centered_OX_data),]
ccc_OY_data <- centered_OY_data[complete.cases(centered_OY_data),]
#print(head(ccc_OX_data))

max_dim = ceiling(length(unique(group$Sample))/2 -1)
print(max_dim)
set.seed(1+1+2020)
if(max_dim<3){
fit_cv_adj <- crossval_o2m_adjR2(X = t(ccc_OX_data), Y = t(ccc_OY_data),a = 2,ax = 0,ay = 0,nr_folds = 2,nr_cores = 4,stripped = TRUE,p_thresh =3000,q_thresh=3000,tol =1e-10,max_iterations =100)
}else{
fit_cv_adj <- crossval_o2m_adjR2(X = t(ccc_OX_data), Y = t(ccc_OY_data),a = 2:max_dim,ax = 0:3,ay = 0:3,nr_folds = 2,nr_cores = 4,stripped =     TRUE,p_thresh =3000,q_thresh=3000,tol =1e-10,max_iterations =100)}

n_n = fit_cv_adj[order(fit_cv_adj$MSE),][1,2]
nx_n =fit_cv_adj[order(fit_cv_adj$MSE),][1,3]
ny_n =fit_cv_adj[order(fit_cv_adj$MSE),][1,4]
fit_o2m <- o2m(X = t(ccc_OX_data), Y = t(ccc_OY_data),n = n_n,nx = nx_n,ny = ny_n)

#t--表示组学X，u--表示组学Y：
#p--表示组学X，q--表示组学Y；
fit_summary = summary(fit_o2m)
fit_R2XY_result <- data.frame(R2X = round(fit_summary$R2_X,4),R2Y = round(fit_summary$R2_Y,4),R2Xjoint = round(fit_summary$R2_Xjoint,4),R2Yjoint =round(fit_summary$R2_Yjoint,4),R2Xhat = round(fit_summary$R2_Xhat,4),R2Yhat = round(fit_summary$R2_Yhat,4),R2Xpred = round(fit_summary$R2_Xpred,4),R2Ypred = round(fit_summary$R2_Ypred,4))
write.table(fit_R2XY_result,"${o}.fit_summary.xls",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
ld_comp <- vector()
sc_comp <- vector()
for(i in seq(1,n_n)){
ld_comp[i] <- paste("pq[",i,"]",sep="")
sc_comp[i] <- paste("tu[",i,"]",sep="")
}
 
X_scores = as.data.frame(fit_o2m$Tt)
#X_scores = as.data.frame(apply(X_scores_raw,2,ScaleUnit))

colnames(X_scores) <- sc_comp
X_scores$Omics = rep(O2O_sp[1],nrow(X_scores))
X_scores$Sample <-row.names(X_scores)
X_scores_new <- merge(X_scores,group,by="Sample",all.x = TRUE,no.dups = FALSE)
rownames(X_scores_new) <- paste(X_scores_new$Sample,"t",sep="_")
Y_scores = as.data.frame(fit_o2m$U)
#Y_scores = as.data.frame(apply(Y_scores_raw,2,ScaleUnit))
#print(head(Y_scores))
colnames(Y_scores) <- sc_comp
Y_scores$Omics = rep(O2O_sp[2],nrow(Y_scores))
Y_scores$Sample <-row.names(Y_scores)
Y_scores_new <- merge(Y_scores,group,by="Sample",all.x = TRUE,no.dups = FALSE)
rownames(Y_scores_new) <- paste(Y_scores_new$Sample,"u",sep="_")
joint_scores <- rbind(X_scores_new,Y_scores_new)

X_loadings = as.data.frame(fit_o2m$W.)
colnames(X_loadings) <- ld_comp
X_loadings$Omics = rep(O2O_sp[1],nrow(X_loadings))
Y_loadings = as.data.frame(fit_o2m$C.)
colnames(Y_loadings) <- ld_comp
Y_loadings$Omics = rep(O2O_sp[2],nrow(Y_loadings))
joint_loadings <- rbind(X_loadings,Y_loadings)
joint_loadings$Omics <- factor(joint_loadings$Omics,levels= c(O2O_sp[1],O2O_sp[2]))

write.table(joint_loadings,"${o}.O2PLS_Loadings.xls",sep="\t",col.names=NA,row.names=TRUE,quote=FALSE)
write.table(joint_scores,"${o}.O2PLS_Scores.xls",sep="\t",col.names=NA,row.names=TRUE,quote=FALSE)

joint_scores$distance <- sqrt(joint_scores[,'tu[1]']^2 + joint_scores[,'tu[1]']^2)
ScaleFactor_raw <- aggregate(joint_scores[,'distance'],by=list(joint_scores[,'Omics']),sum)
ScaleFactor <- ScaleFactor_raw[ScaleFactor_raw==O2O_sp[1],2]/ScaleFactor_raw[ScaleFactor_raw==O2O_sp[2],2]

joint_scores$`Scaled tu[1]` <- ifelse(joint_scores$Omics==O2O_sp[1],joint_scores[,'tu[1]']*1,joint_scores[,'tu[1]']*ScaleFactor)
joint_scores$`Scaled tu[2]` <- ifelse(joint_scores$Omics==O2O_sp[1],joint_scores[,'tu[2]']*1,joint_scores[,'tu[2]']*ScaleFactor)
"""
f = Template(Rcmd,output_encoding='utf-8',input_encoding='utf-8')
Rcmd = f.render(omics_x = args.omics_x,omics_y = args.omics_y,group = args.group,oxoy = args.oxoy,scale = args.scale,path=args.path ,o=args.out,encoding="utf-8",log_x=args.log_x,log_y=args.log_y)
fout=open('%s.O2PLS.cmd.r' % args.out,'w')
print('ok')
fout.write(Rcmd)
fout.close()




