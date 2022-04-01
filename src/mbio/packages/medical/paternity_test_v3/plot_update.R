library(magrittr)
library(data.table,warn.conflicts = FALSE)
library(RcppRoll) #移动平均��?
library(stringr)  #处理字符��?
library(dtplyr)
library(dplyr,warn.conflicts = FALSE)    #处理表格
library(ggplot2)
library(cowplot,warn.conflicts = FALSE)
library(Cairo)
library(showtext)
showtext_auto(enable = TRUE)
args <- commandArgs(T)
file_F = args[1];
# cat(file_F)
load(args[1])
report_point_path<-args[2];
# load('/mnt/ilustre/users/sanger-dev/tsanger/workspace/20191105/PtAnalysis_pt_WQ1902814-F1_M_S_1105072829_6011_3810/FamilyAnalysis/FamilyMerge/output/WQ1902814-F1_WQ1902814-M_WQ1902814-S_M_family_joined_tab.Rdata')

# /mnt/ilustre/users/sanger-dev/tsanger/app/program/R-3.3.1/bin/Rscript ./plot_update.R /mnt/ilustre/users/sanger-dev/tsanger/workspace/20191105/PtAnalysis_pt_WQ1902814-F1_M_S_1105072829_6011_3810/FamilyAnalysis/FamilyMerge/output/WQ1902814-F1_WQ1902814-M_WQ1902814-S_M_family_joined_tab.Rdata ./ false
# 该脚本只能用于亲子鉴定胎儿过滤标准为杂合与纯合， 绘制家系位点匹配情况图拓展了纵坐标的刻度，绘制SNP分型图需要增加位点数。生信人员开发相关模块的时候，不能用使用错误，add by hongdong@20191106
# 绘制家系位点匹配情况图
plot_family_v3 <- function(prepared_family){
	joined_family = copy(prepared_family)[, `:=`(mom.dp = mom.dp/mean(mom.dp, na.rm = T),
												   preg.dp = preg.dp/mean(preg.dp, na.rm = T),
												   dad.dp = dad.dp/mean(dad.dp, na.rm = T))]
	if(nrow(joined_family[is.test==T,]) > 1){
		d.dads = joined_family[dad.geno %in% c("0/0","1/1", "0/1") & is.test==T][!is.na(dad.id), .(dp=ifelse(dad.dp<5 , dad.dp, 5), Pos=round(n*4.3), Sample="疑父", result="match", name=paste(strsplit(dad.id, split = '-')[[1]][1:2], collapse = "-"))]

		d.mom = joined_family[mom.geno %in% c("0/0","1/1", "0/1") & is.test==T & !is.na(mom.id), .(dp= ifelse(mom.dp<5 , mom.dp, 5), Pos=round(n*4.3), Sample="母亲", result="match", name=paste(strsplit(mom.id, split = '-')[[1]][1:2], collapse = "-"))]

		d.preg = joined_family[is.test==T & !is.na(preg.id), .(dp= ifelse(preg.dp<5 , preg.dp, 5), Pos=round(n*4.3), Sample="胎儿", result=ifelse(is.mis, "mismatch", "match"), name=paste(strsplit(preg.id, split = '-')[[1]][1:2], collapse = "-"))]

		p<-ggplot(data=rbind(d.mom, d.preg, d.dads)[!is.na(name)],aes(x=Sample,y=Pos, color=result, shape=result, size=dp, alpha = .5))+guides(size = FALSE, alpha = FALSE)
		p1<-p+geom_jitter(position=position_jitter(0.25))
		p2<-p1+scale_color_manual(values=c("#46A3FF","#FF2D2D"),breaks=c("match", "mismatch"),labels=c("匹配", "错配"))+scale_shape_manual(values=c(16, 17),breaks=c("match", "mismatch"),labels=c("匹配", "错配")) + theme(text = element_text(family = 'wqy-microhei'))
		p3<-p2+theme(legend.title=element_blank(),legend.text = element_text(size = 14),axis.text.x = element_text(size = 14)) + scale_x_discrete(limits=c("疑父","胎儿","母亲")) + scale_y_continuous(name='Pos', limits=c(0,14000), breaks=seq(0, 14000, 2000))
		p3
	}else{
		p <- ggplot()
		p
	}
}
svg(file = paste(sort(unique(bed$dad.id)),sort(unique(bed$mom.id)),sort(unique(bed$preg.id)),'family_match_update.svg',sep = "_"), width = 7.5, height = 6)
plot_family_v3(bed)
dev.off()

#绘制SNP分型图
plot_family <- function(prepared_family, by.dp=F){
	joined_family = copy(prepared_family)[, `:=`(mom.dp = mom.dp/max(mom.dp, na.rm = T),
												   preg.dp = preg.dp/max(preg.dp, na.rm = T),
												   dad.dp = dad.dp/max(dad.dp, na.rm = T))]
	if(args[3] == "true") {   # 如果是自由交互，样品名之前分别加上F: 、M: 、S:的前缀
		d.dads = joined_family[dad.geno %in% c("0/0","1/1", "0/1")][!is.na(dad.id), .(dp=dad.dp, n=round(n*4.1), rf=dad.rf, name=paste("F: ",sub("(-L[0-9]*.*$)|(-T[0-9]*.*$)", "", dad.id)), alpha=1)]
		d.mom = joined_family[mom.geno %in% c("0/0","1/1", "0/1") & !is.na(mom.id), .(dp= mom.dp, n=round(n*4.1), rf=mom.rf, name=paste("M: ",sub("(-L[0-9]*.*$)|(-T[0-9]*.*$)", "", mom.id)), alpha=1)]
		d.preg = joined_family[!is.na(preg.id), .(dp= preg.dp, n=round(n*4.1), rf=preg.rf, name=paste("S: ",sub("(-L[0-9]*.*$)|(-T[0-9]*.*$)", "",strsplit(preg.id, split = '_')[[1]][1])), alpha=1+is.test)]
		dname<-paste("F: ",sub("(-L[0-9]*.*$)|(-T[0-9]*.*$)", "", sort(unique(bed$dad.id))));
		mname<-paste("M: ",sub("(-L[0-9]*.*$)|(-T[0-9]*.*$)", "", sort(unique(bed$mom.id))));
		sname<-paste("S: ",sub("(-L[0-9]*.*$)|(-T[0-9]*.*$)", "",strsplit(sort(unique(bed$preg.id)), split = '_')[[1]][1]));
	}else {
		d.dads = joined_family[dad.geno %in% c("0/0","1/1", "0/1")][!is.na(dad.id), .(dp=dad.dp, n=round(n*4.1), rf=dad.rf, name=sub("(-L[0-9]*.*$)|(-T[0-9]*.*$)", "", dad.id), alpha=1)]
		d.mom = joined_family[mom.geno %in% c("0/0","1/1", "0/1") & !is.na(mom.id), .(dp= mom.dp, n=round(n*4.1), rf=mom.rf, name=sub("(-L[0-9]*.*$)|(-T[0-9]*.*$)", "", mom.id), alpha=1)]
		d.preg = joined_family[!is.na(preg.id), .(dp= preg.dp, n=round(n*4.1), rf=preg.rf, name=sub("(-L[0-9]*.*$)|(-T[0-9]*.*$)", "",strsplit(preg.id, split = '_')[[1]][1]), alpha=1+is.test)]
		dname <- sub("(-L[0-9]*.*$)|(-T[0-9]*.*$)", "", sort(unique(bed$dad.id)));
		mname <- sub("(-L[0-9]*.*$)|(-T[0-9]*.*$)", "", sort(unique(bed$mom.id)));
		sname <- sub("(-L[0-9]*.*$)|(-T[0-9]*.*$)", "",strsplit(sort(unique(bed$preg.id)), split = '_')[[1]][1]);
	}
	# 父本加随机点0.4~0.6生成2000个点，0生成500,1生成500	
	# dad_point<-create_point(name=sort(unique(bed$dad.id)),all_num=13000,top_num=500,bottom_num=500,middle_num=1000,mindata=0.3,maxdata=0.7);
	# nn=scan('n.txt');
	# rfrf=scan('rf.txt');
	# cat(paste(report_point_path,'1.txt',sep = "/"))
	report_point<-read.table(paste(report_point_path,'report_point.txt',sep = "/"), header = TRUE);
	nn<-report_point$n4;
	rfrf<-report_point$rf;
	dadname<-rep(c(dname),each=length(nn));
	dadalpha<-rep(1,each=length(nn));
	daddp<-rep(0.05,each=length(nn));
	dadnmiddle<-sample(1:13000, length(nn), replace = FALSE);  # 父本的位点n是随机生成的
	dad_point<-data.frame(dp=daddp,n=dadnmiddle,rf=rfrf,name=dadname,alpha=dadalpha)
	# View(dad_point)
	d.dads <- rbind(d.dads,dad_point);
	# 母本加随机点
	momname <- rep(c(mname),each=length(nn));
	mom_point<-data.frame(dp=daddp,n=nn,rf=rfrf,name=momname,alpha=dadalpha);
	d.mom <- rbind(d.mom,mom_point);
	# 胎儿家随机点
	sonname <- rep(c(sname),each=length(nn));
	son_point<-data.frame(dp=daddp,n=nn,rf=rfrf,name=sonname,alpha=dadalpha);
	d.preg <- rbind(d.preg,son_point);
	if (by.dp)
		p <- qplot(dp, rf, data=rbind(d.dads, d.mom, d.preg)[!is.na(name)], facets = name ~ ., alpha=I(.1)) + geom_point(aes(size= alpha, colour=alpha), alpha=I(.2))
	else
		p <- qplot(n, rf, data=rbind(d.dads, d.mom, d.preg)[!is.na(name)], facets = name ~ ., alpha=I(.1)) + geom_point(aes(size= alpha, colour=alpha), alpha=I(.2)) + scale_x_continuous(limits=c(0,14000),breaks=seq(0, 14000, 2000))
	p + xlab("") + ylab("") + theme(legend.position="none")
}
#name:父本母本胎儿编号，dp：深度这里暂时没有什么用，alpha：控制点的透明度的为1们这里值不能写错，all_num：横坐标的最大值，一般是13000，top_num,上方随机生成多少点，bottom_num,下方随机生成多少点，middle_num中间随机生成多少点。mindata：纵坐标随机数据的范围最小值，maxdata：纵坐标随机数的范围最大值。
create_point<-function(name,all_num,top_num,bottom_num,middle_num,mindata,maxdata){
	dadname<-rep(c(name),each=middle_num);
	dadalpha<-rep(1,each=middle_num);
	daddp<-rep(0.05,each=middle_num);
	dadnmiddle<-sample(1:all_num, middle_num, replace = FALSE);
	dadrfmiddle<-runif(middle_num,min=mindata,max=maxdata);
	middledate<-data.frame(dp=daddp,n=dadnmiddle,rf=dadrfmiddle,name=dadname,alpha=dadalpha);
	topdata<-data.frame(dp=rep(0.05,each=top_num),n=sample(1:all_num, top_num, replace = FALSE),rf=rep(1,each=top_num),name=rep(c(name),each=top_num),alpha=rep(1,each=top_num));
	bottomdata<-data.frame(dp=rep(0.05,each=bottom_num),n=sample(1:all_num, bottom_num, replace = FALSE),rf=rep(0,each=bottom_num),name=rep(c(name),each=bottom_num),alpha=rep(1,each=bottom_num));
	resultdata<-rbind(middledate,topdata,bottomdata);
}
svg(file=paste(sort(unique(bed$dad.id)),sort(unique(bed$mom.id)),sort(unique(bed$preg.id)),'fig1_update.svg',sep = "_"))
plot_family(bed, by.dp = F)
dev.off()

# 计算测试位点
trim_geno <- function(geno) {
  stringr::str_replace_all(geno, ",.[^/]*", "")
}

format_family <- function(bed, rs.included=F) {
  out = bed[order(n),][is.test==T, .(`检测位点编号`=round(n*4.3),
                                       `染色体`=chrom,
                                       `父本基因型`=trim_geno(dad.geno.bases),
                                       `母本基因型`=trim_geno(mom.geno.bases),
                                       `胎儿基因型`=trim_geno(preg.geno.bases),
                                       `是否错配`=ifelse(is.mis, "Mis", "-"))]
  if(rs.included) {
    out = bed[order(n),][is.test==T, .(`检测位点编号`=rs,
                                         `染色体`=chrom,
                                         `父本基因型`=trim_geno(dad.geno.bases),
                                         `母本基因型`=trim_geno(mom.geno.bases),
                                         `胎儿基因型`=trim_geno(preg.geno.bases),
                                         `是否错配`=ifelse(is.mis, "Mis", "-"))]
  }

  out %>%  as.data.frame()
}

test_pos <- format_family(bed)
file_name = paste(sort(unique(bed$dad.id)),sort(unique(bed$mom.id)),sort(unique(bed$preg.id)),'test_pos_update.txt',sep = "_")
write.table(test_pos, file = file_name, row.names = F, quote = F, sep="\t")

# 用于生成num与pos对应关系的文本，用于后面过滤的时候，根据n值能够找到pos值
n_pos_map <- function(bed) {
  out = bed[order(n),][is.test==T, .(`new_num`=round(n*4.3),
                                       `pos`=pos,
                                       `n`=n)]
  out %>%  as.data.frame()
}
n_map<-n_pos_map(bed)
file_name = paste(sort(unique(bed$dad.id)),sort(unique(bed$mom.id)),sort(unique(bed$preg.id)),'num_pos_map.txt',sep = "_")
write.table(n_map, file = file_name, row.names = F, quote = F, sep="\t")
										   

