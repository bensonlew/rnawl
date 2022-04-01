#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/5/23 12:49
@file    : dia_quality_control.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""


import sys
import os
import pandas as pd
from collections import OrderedDict
import argparse

reload(sys)
sys.setdefaultencoding("utf-8")

parser = argparse.ArgumentParser(description="dia流程的质控处理")
parser.add_argument("-report", type=str, required=True, help="dia流程结果里的report的xls文件，一般就特别大的那个")
parser.add_argument("-exp_ana_file", type=str, required=True, help="dia结果的*ExperimentAnalysis.txt")
parser.add_argument("-protein", type=str, required=True, help='蛋白表，跟PD搜库一样的表格')
parser.add_argument("-normal", type=str, required=True, help='包含蛋白coverage信息的文件，名字一般为：Report_Protein quant Cov NSSI (Normal).xls')
args = parser.parse_args()
# if len(sys.argv) != 3:
#     exit('%s report.xls exp_ana_file'%sys.argv[0])

# report = sys.argv[1]
report = args.report
exp_ana_file = args.exp_ana_file
protein_file = args.protein
normal_file = args.normal
try:
    report_df = pd.read_csv(report, sep='\t')
except:
    report_df = pd.read_excel(report, sep='\t')

# plot Peptide length distribution
lw = open('Peptide_sequence_length.list', 'w')
peps = list()
for pep in report_df['PEP.GroupingKey']:
    pep = ''.join(list(filter(str.isalpha, pep)))
    if pep in peps:
        continue
    if len(pep):
        lw.write(pep+'\t'+str(len(pep))+'\n')
        peps.append(pep)
lw.close()
r_cmd = r'''options(warn=-100)

his <- read.table("Peptide_sequence_length.list", header = T, sep = "\t")
dat <- as.data.frame(ftable(his[,2]))
colnames(dat) <- c("Length", "Number")

write.table(dat, file = "Peptide_length_distribution.xls", quote = F, sep = "\t", col.names = T, row.names = F)

pdf("Peptide_length_distribution.pdf", w=10, h=7)
par(mar=c(5.1,6.1,4.1,2.1))
bp <- barplot(dat$Number, col = "#4F81BD", border=F, las=0,width=1,space=0.5, plot=T,beside=T, axisnames=F, xlab = "Peptide length", ylab = "Peptide number", main = "Peptide length distribution", cex.main=1.2,cex.axis=0.8,las=1)
#bp <- barplot(dat$Number, col = "#4F81BD", border=F, las=0,width=1,space=0.5,ylim=c(-0.02*max(dat$Number),max(dat$Number)*1.2), plot=T,beside=T, axisnames=F, xlab = "Peptide length", ylab = "Peptide number", main = "Peptide length distribution", cex.main=1.2,cex.axis=0.8,las=1)

text(x=bp,y=-0.013*max(dat$Number),srt=70,adj=1,labels=dat$Length,xpd=T,cex=0.8)

jump<-max(dat$Number)/10
text(bp,(dat$Number+jump/3),labels=dat$Number,xpd=T,cex=0.5)
dev.off()
'''
with open('Peptide_sequence_length.r', 'w') as rw:
    rw.write(r_cmd)
os.system('Rscript Peptide_sequence_length.r')


# plot Peptide number distribution
pro2peps = dict()
for pro, pep in zip(report_df['PG.ProteinGroups'], report_df['PEP.GroupingKey']):
    pros = pro.split(';')
    pep = ''.join(list(filter(str.isalpha, pep)))
    for p in pros:
        if not p in pro2peps:
            pro2peps[p] = list()
        if pep not in pro2peps[p]:
            pro2peps[p].append(pep)
with open('Peptide_number.txt', 'w') as pnw:
    for pro in pro2peps:
        l = len(pro2peps[pro])
        if l:
            pnw.write(pro+'\t'+str(l)+'\n')
r_cmd = r'''options(warn=-100)

df_raw <- read.table("Peptide_number.txt", header = T, sep = "\t")
df <- df_raw[,2]
for (i in 1:10){
	if (length(which(df <10*i)) >= 0.9*length(df)){ 
		df1<-length(which(df>0*i&df<=1*i))
		df2<-length(which(df>1*i&df<=2*i))
		df3<-length(which(df>2*i&df<=3*i))
		df4<-length(which(df>3*i&df<=4*i))
		df5<-length(which(df>4*i&df<=5*i))
		df6<-length(which(df>5*i&df<=6*i))
		df7<-length(which(df>6*i&df<=7*i))
		df8<-length(which(df>7*i&df<=8*i))
		df9<-length(which(df>8*i&df<=9*i))
		df10<-length(which(df>9*i&df<=10*i))
		df11<-length(which(df>10*i&df<=11*i))
		df12<-length(which(df>11*i&df<=12*i))
		df13<-length(which(df>12*i&df<=13*i))
		df14<-length(which(df>13*i&df<=14*i))
		df15<-length(which(df>14*i&df<=15*i))
		df16<-length(which(df>15*i))

		num <- c(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16)
		dis <- as.character(c(paste(0*i+1,'-',1*i),paste(1*i+1,'-',2*i),paste(2*i+1,'-',3*i),paste(3*i+1,'-',4*i),paste(4*i+1,'-',5*i),paste(5*i+1,'-',6*i),paste(6*i+1,'-',7*i),paste(7*i+1,'-',8*i),paste(8*i+1,'-',9*i),paste(9*i+1,'-',10*i),paste(10*i+1,'-',11*i),paste(11*i+1,'-',12*i),paste(12*i+1,'-',13*i),paste(13*i+1,'-',14*i),paste(14*i+1,'-',15*i),paste('>',15*i)))
		if (i == 1) {
			dis <- as.character(c(1:15, ">15"))
		}
		dat <- data.frame(dis,num)
		colnames(dat) <- c("Peptide_number", "Protein_number")
		write.table(dat, file = "Peptide_number_distribution.xls", quote = F, sep = "\t", col.names = T, row.names = F)

		pdf("Peptide_number_distribution.pdf", w=10,h=7)
		par(mar=c(5.1,6.1,4.1,2.1))
		bp <- barplot(num, col="#4F81BD",border=FALSE, axisnames=F,las=0,width=1,space=0.5, plot=T,beside=TRUE,main="Peptide number distribution",cex.main=1.2,xlab="Peptide number",ylab="Protein number",cex.axis=0.8,las=1)
		#bp <- barplot(num, col="#4F81BD",border=FALSE,ylim=c(-0.02*max(num),max(num)*1.2), axisnames=F,las=0,width=1,space=0.5, plot=T,beside=TRUE,main="Peptide number distribution",cex.main=1.2,xlab="Peptide number",ylab="Protein number",cex.axis=0.8,las=1)

		text(x=bp,y=-0.013*max(num),srt=70,adj=1,labels=dis,xpd=T,cex=0.8)

		jump<-max(num)/10
		text(bp,(num+jump/3),labels=num,xpd=T,cex=0.8)

		dev.off()
		break
	}
}
'''
with open('Peptide_number_distribution.r', 'w') as rw:
    rw.write(r_cmd)
os.system('Rscript Peptide_number_distribution.r')

# plot dMass
dw = open('dMass.txt', 'w')
for a,b in zip(report_df['FG.PrecMz'], report_df['FG.PrecMzCalibrated']):
    if a and b:
        dw.write(str(b)+'\t'+str((b-a)/b*1000000)+'\n')
dw.close()
r_cmd = r'''dat<-read.table("dMass.txt",sep="\t",header=T)
df<-data.frame(dat[,1], dat[,2])
pdf("dMass.pdf", 10, 8)
par(mar=c(4.1, 5.1, 4.1, 2))
# plot(df, pch=".", col = "red", bty="o", ylim = c(-20,20), xlab = "m/z(Da)", ylab = "DeltaM(ppm)")
plot(df, pch=".", col = "red", bty="o", ylim = c(0,5), xlab = "m/z(Da)", ylab = "DeltaM(ppm)")
dev.off()
'''
with open('dMass.r', 'w') as rw:
    rw.write(r_cmd)
os.system('Rscript dMass.r')
del report_df


# plot protein information
# exp_ana_file = sys.argv[2]
info2num = OrderedDict()
with open(exp_ana_file, 'r') as er:
    for line in er:
        if not line.strip():
            continue
        if u'Sparse Profiles - ' in line:
            info, num = line.strip().split('\t')
            info = info.split('Sparse Profiles - ')[1]
            num = num.split(' of ')[0]
            info2num[info] = num
info_df = pd.DataFrame([info2num])
tmp = info_df.columns.tolist()
tmp[-1], tmp[-2] = tmp[-2], tmp[-1]
info_df = info_df[tmp]
info_df.to_csv('Protein_information.xls', sep='\t', header=True, index=False)
r_cmd = r'''options(warn=-100)

df_raw <- read.table("Protein_information.xls", header = T, sep = "\t", check.names = F)
df <- as.data.frame(t(df_raw))
names <- rownames(df)
number <- df$V1

dat <- as.data.frame(names, number)

pdf("Protein_information.pdf", w=13,h=9)
par(mar=c(5.1,6.1,4.1,2.1))

bp <- barplot(number, col="#4F81BD", names.arg = names, border=FALSE,las=0,width=1,space=0.5,ylim=c(-0.02*max(number),max(number)*1.2), plot=T,beside=TRUE,main="Protein information",cex.main=1.2,xlab="",ylab="Number",cex.axis=0.8,las=1)

jump<-max(number)/10
text(bp,(number+jump/3),labels=number,xpd=T,cex=0.8)

dev.off()
'''
with open('Protein_information.r', 'w') as rw:
    rw.write(r_cmd)
os.system('Rscript Protein_information.r')

# plot Protein molecular weight distribution
try:
    protein_df = pd.read_csv(protein_file, sep='\t')
except:
    protein_df = pd.read_excel(protein_file, sep='\t')
# pro2weight = dict()
# for pro, weight in zip(report_df['PG.ProteinGroups'], report_df['FG.Quantity']):
#     pros = pro.split(';')
#     weight = weight / 1000
#     for p in pros:
#         if not p in pro2weight:
#             pro2weight[p] = list()
#         pro2weight[p].append(weight)
# with open('Protein_MW.txt', 'w') as pmw:
#     for pro in pro2weight:
#         weight = pro2weight[pro]
#         if not weight:
#             continue
#         weight = sum(weight) / len(weight)
#         pmw.write(pro + '\t' + str(weight) + '\n')
weight_df = protein_df[['Accession', 'MW [kDa]']]
weight_df.to_csv('Protein_MW.txt', sep='\t', header=True, index=False)
r_cmd = r'''options(warn=-100)

df_raw <- read.table("Protein_MW.txt", header = T, sep = "\t")
df <- df_raw[,2]
rank = 10
for (i in 1:10){
	if (length(which(df <100*i)) >= 0.9*length(df)){ 
		df1<-length(which(df>=0&df<=1*rank*i))
		df2<-length(which(df>1*rank*i&df<=2*rank*i))
		df3<-length(which(df>2*rank*i&df<=3*rank*i))
		df4<-length(which(df>3*rank*i&df<=4*rank*i))
		df5<-length(which(df>4*rank*i&df<=5*rank*i))
		df6<-length(which(df>5*rank*i&df<=6*rank*i))
		df7<-length(which(df>6*rank*i&df<=7*rank*i))
		df8<-length(which(df>7*rank*i&df<=8*rank*i))
		df9<-length(which(df>8*rank*i&df<=9*rank*i))
		df10<-length(which(df>9*rank*i&df<=10*rank*i))
		df11<-length(which(df>10*rank*i&df<=11*rank*i))
		df12<-length(which(df>11*rank*i&df<=12*rank*i))
		df13<-length(which(df>12*rank*i&df<=13*rank*i))
		df14<-length(which(df>13*rank*i&df<=14*rank*i))
		df15<-length(which(df>14*rank*i&df<=15*rank*i))
		df16<-length(which(df>15*rank*i))

		mw <- c(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16)
		dis <- as.character(c(paste(0*rank*i+1,'-',1*rank*i),paste(1*rank*i+1,'-',2*rank*i),paste(2*rank*i+1,'-',3*rank*i),paste(3*rank*i+1,'-',4*rank*i),paste(4*rank*i+1,'-',5*rank*i),paste(5*rank*i+1,'-',6*rank*i),paste(6*rank*i+1,'-',7*rank*i),paste(7*rank*i+1,'-',8*rank*i),paste(8*rank*i+1,'-',9*rank*i),paste(9*rank*i+1,'-',10*rank*i),paste(10*rank*i+1,'-',11*rank*i),paste(11*rank*i+1,'-',12*rank*i),paste(12*rank*i+1,'-',13*rank*i),paste(13*rank*i+1,'-',14*rank*i),paste(14*rank*i+1,'-',15*rank*i),paste('>',15*rank*i)))
		dat <- data.frame(dis,mw)
		colnames(dat) <- c("Molecular_weight(kDa)", "Protein_number")
		write.table(dat, file = "Protein_molecular_weight_distribution.xls", quote = F, sep = "\t", col.names = T, row.names = F)

		pdf("Protein_molecular_weight_distribution.pdf", w=10,h=7)
		par(mar=c(5.1,6.1,4.1,2.1))
		bp <- barplot(mw, col="#4F81BD",border=FALSE, axisnames=F,las=0,width=1,space=0.5,ylim=c(-0.02*max(mw),max(mw)*1.2), plot=T,beside=TRUE,main="Protein molecular weight distribution",cex.main=1.2,xlab="Protein molecular weight (kDa)",ylab="Protein number",cex.axis=0.8,las=1)

		text(x=bp,y=-0.013*max(mw),srt=70,adj=1,labels=dis,xpd=T,cex=0.8)

		jump<-max(mw)/10
		text(bp,(mw+jump/3),labels=mw,xpd=T,cex=0.8)

		dev.off()
		break
	}
}
'''
with open('Protein_molecular_weight_distribution.r', 'w') as rw:
    rw.write(r_cmd)
os.system('Rscript Protein_molecular_weight_distribution.r')

# plot Protein Coverage distribution

# coverage_df = protein_df[['Accession', 'Coverage [%]']]
# coverage_df.to_csv('Protein_Coverage.txt', sep='\t', header=True, index=False)

try:
    normal_df = pd.read_csv(normal_file, sep='\t')
except:
    normal_df = pd.read_excel(normal_file, sep='\t')
pro2cover = dict()
for p,c in zip(normal_df['PG.ProteinAccessions'], normal_df['PG.Coverage']):
    p = p.split(';')
    c = [float(x.strip('%')) for x in c.split(';')]
    for p_,c_ in zip(p,c):
        if p_ not in pro2cover:
            pro2cover[p_] = list()
        pro2cover[p_].append(c_)
with open('Protein_Coverage.txt', 'w') as cw:
    cw.write('Accession\tCoverage [%]\n')
    for p in pro2cover:
        cs = pro2cover[p]
        if cs:
            cw.write(p+'\t'+str(max(cs))+'\n')

r_cmd = '''options(warn=-100)
library("RColorBrewer")

df_raw <- read.table("Protein_Coverage.txt", header = T, sep = "\t")
df <- df_raw[,2]
#df1<-length(which(df>0&df<=1)) 20180129 yanjun
df1<-length(which(df<=1)) 
df2<-length(which(df>1&df<=5))
df3<-length(which(df>5&df<=10))
df4<-length(which(df>10&df<=20))
df5<-length(which(df>20&df<=40))
df6<-length(which(df>40&df<=60))
df7<-length(which(df>60&df<=80))
df8<-length(which(df>80&df<=100))

cov <- c(df1,df2,df3,df4,df5,df6,df7,df8)
dis <- as.character(c("<1","1-5","6-10","11-20","21-40","41-60","61-80",">80"))
dat <- data.frame(dis,cov)
colnames(dat) <- c("Coverage(%)", "Protein_number")
write.table(dat, file = "Protein_coverage_distribution.xls", quote = F, sep = "\t", col.names = T, row.names = F)

pdf("Protein_coverage_distribution.pdf", w=10, h=8)
par(mar=c(5.1,6.1,4.1,2.1))
pie(cov, col=brewer.pal(8,"Set3"), label=cov, clockwise=T, cex=0.8, radius=0.85, border="white")
legend("right", dis, title="Coverage(%)", fill=brewer.pal(8,"Set3"), border=F, bty="n")

dev.off()
'''
with open('Protein_coverage_distribution.r', 'w') as rw:
    rw.write(r_cmd)
os.system('Rscript Protein_coverage_distribution.r')
print('all done')