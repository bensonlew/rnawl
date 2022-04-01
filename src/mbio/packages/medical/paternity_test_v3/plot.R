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
# /mnt/ilustre/users/sanger-dev/app/program/R-3.3.1/bin/Rscript plot.R /mnt/ilustre/users/sanger-dev/workspace/20161216/PtProcess_pt_process_test/output/family_joined_tab.Rdata
# 修改了父本 母本 胎儿的名字
# load('D:/medical/code_pt/tool/family_joined_tab.Rdata') #加载bed
# load('family_joined_tab.Rdata') #加载bed

args <- commandArgs()
file_F = args[6]
load(args[6])
###补充胎儿浓度的计��?
baby = bed[mom.geno=="1/1" & str_length(mom.geno.bases)==3 & preg.dp >= 10,][preg.dp >= quantile(preg.dp, 0.1) & preg.rf <= quantile(preg.rf, 0.99)]

d0 = baby[,.(ref=preg.ref.dp, all=preg.dp)]

m = function(para)
{
  d.error = dbinom(x=d0$ref, size = d0$all, prob = para[2])
  d.preg = dbinom(x=d0$ref, size = d0$all, prob = para[3])
  d.all=para[1]*d.error+(1-para[1])*d.preg
  max.like=sum(log(d.all))
  return(-max.like)
}

#非线性最小化的函数nlminb
model.est = nlminb(c(.9, .003, .1),
                   m,
                   lower = c(0.001, 0.0001, 0.01),
                   upper = c(0.999, 0.005, 0.3))

baby_per = c(p.baby = model.est$par[3], p.error = model.est$par[2], pp = 1-model.est$par[1])

#胎儿浓度
percent = (baby_per[1]*2*100) %>% round(2)
#测序错误
error = (baby_per[2]*100)%>% round(4)
#胎儿信号占比
s_signal = (baby_per[3]*100)%>% round(4)
#胎儿深度
dp_preg = mean(bed$preg.dp,na.rm=T)%>%round(2)
#母本深度
dp_mom = mean(bed$mom.dp,na.rm=T) %>% round(2)
#母本与胎儿匹配程度
mom_preg = bed[mom.geno=="0/0" & preg.dp > 20, 100* length(preg.rf[preg.rf > 0.7])/length(preg.rf)] %>% round(4)


plot_family <- function(prepared_family, by.dp=F){
  joined_family = copy(prepared_family)[, `:=`(mom.dp = mom.dp/max(mom.dp, na.rm = T),
                                               preg.dp = preg.dp/max(preg.dp, na.rm = T),
                                               dad.dp = dad.dp/max(dad.dp, na.rm = T))]
  if(args[7] == "true") {   # 如果是自由交互，样品名之前分别加上F: 、M: 、S:的前缀
    d.dads = joined_family[dad.geno %in% c("0/0","1/1", "0/1")][!is.na(dad.id), .(dp=dad.dp, n, rf=dad.rf, name=paste("F: ",sub("(-L[0-9]*.*$)|(-T[0-9]*.*$)", "", dad.id)), alpha=1)]
  d.mom = joined_family[mom.geno %in% c("0/0","1/1", "0/1") & !is.na(mom.id), .(dp= mom.dp, n, rf=mom.rf, name=paste("M: ",sub("(-L[0-9]*.*$)|(-T[0-9]*.*$)", "", mom.id)), alpha=1)]
  d.preg = joined_family[!is.na(preg.id), .(dp= preg.dp, n, rf=preg.rf, name=paste("S: ",sub("(-L[0-9]*.*$)|(-T[0-9]*.*$)", "",strsplit(preg.id, split = '_')[[1]][1])), alpha=1+is.test)]
  }else {
  d.dads = joined_family[dad.geno %in% c("0/0","1/1", "0/1")][!is.na(dad.id), .(dp=dad.dp, n, rf=dad.rf, name=sub("(-L[0-9]*.*$)|(-T[0-9]*.*$)", "", dad.id), alpha=1)]
  #View(d.dads)
  d.mom = joined_family[mom.geno %in% c("0/0","1/1", "0/1") & !is.na(mom.id), .(dp= mom.dp, n, rf=mom.rf, name=sub("(-L[0-9]*.*$)|(-T[0-9]*.*$)", "", mom.id), alpha=1)]
  #View(d.mom)
  d.preg = joined_family[!is.na(preg.id), .(dp= preg.dp, n, rf=preg.rf, name=sub("(-L[0-9]*.*$)|(-T[0-9]*.*$)", "",strsplit(preg.id, split = '_')[[1]][1]), alpha=1+is.test)]
  #View(d.preg)
  }

  if (by.dp)
    p <- qplot(dp, rf, data=rbind(d.dads, d.mom, d.preg)[!is.na(name)], facets = name ~ ., alpha=I(.1)) + geom_point(aes(size= alpha, colour=alpha), alpha=I(.2))
  else
    p <- qplot(n, rf, data=rbind(d.dads, d.mom, d.preg)[!is.na(name)], facets = name ~ ., alpha=I(.1)) + geom_point(aes(size= alpha, colour=alpha), alpha=I(.2))
  p + xlab("") + ylab("") + theme(legend.position="none")
}

# 家系��?
# svg(file = paste(bed$dad.id,bed$mom.id,bed$preg.id,'family.svg',sep = "_"))
svg(file = paste(sort(unique(bed$dad.id)),sort(unique(bed$mom.id)),sort(unique(bed$preg.id)),'family.svg',sep = "_"))
plot_family(bed, by.dp = T)
dev.off()

#图1
# svg(file=paste(bed$dad.id,bed$mom.id,bed$preg.id,'fig1.svg',sep = "_"))
svg(file=paste(sort(unique(bed$dad.id)),sort(unique(bed$mom.id)),sort(unique(bed$preg.id)),'fig1.svg',sep = "_"))
plot_family(bed, by.dp = F)
dev.off()

# 亲子V3新增家系位点匹配图（V3）
plot_family_v3 <- function(prepared_family){
  joined_family = copy(prepared_family)[, `:=`(mom.dp = mom.dp/mean(mom.dp, na.rm = T),
                                               preg.dp = preg.dp/mean(preg.dp, na.rm = T),
                                               dad.dp = dad.dp/mean(dad.dp, na.rm = T))]
  if(nrow(joined_family[is.test==T,]) > 1){
  d.dads = joined_family[dad.geno %in% c("0/0","1/1", "0/1") & is.test==T][!is.na(dad.id), .(dp=ifelse(dad.dp<5 , dad.dp, 5), Pos=n, Sample="疑父", result="match", name=paste(strsplit(dad.id, split = '-')[[1]][1:2], collapse = "-"))]

  d.mom = joined_family[mom.geno %in% c("0/0","1/1", "0/1") & is.test==T & !is.na(mom.id), .(dp= ifelse(mom.dp<5 , mom.dp, 5), Pos=n, Sample="母亲", result="match", name=paste(strsplit(mom.id, split = '-')[[1]][1:2], collapse = "-"))]

  d.preg = joined_family[is.test==T & !is.na(preg.id), .(dp= ifelse(preg.dp<5 , preg.dp, 5), Pos=n, Sample="胎儿", result=ifelse(is.mis, "mismatch", "match"), name=paste(strsplit(preg.id, split = '-')[[1]][1:2], collapse = "-"))]

  p<-ggplot(data=rbind(d.mom, d.preg, d.dads)[!is.na(name)],aes(x=Sample,y=Pos, color=result, shape=result, size=dp, alpha = .5))+guides(size = FALSE, alpha = FALSE)
  p1<-p+geom_jitter(position=position_jitter(0.25))
  p2<-p1+scale_color_manual(values=c("#32CD32","#FF6666"),breaks=c("match", "mismatch"),labels=c("匹配", "错配"))+scale_shape_manual(values=c(16, 17),breaks=c("match", "mismatch"),labels=c("匹配", "错配")) + theme(text = element_text(family = 'wqy-microhei'))
  p3<-p2+theme(legend.title=element_blank(),legend.text = element_text(size = 14),axis.text.x = element_text(size = 14)) + scale_x_discrete(limits=c("疑父","胎儿","母亲"))
  p3}else{
  	p <- ggplot()
  	p
  }

}

svg(file = paste(sort(unique(bed$dad.id)),sort(unique(bed$mom.id)),sort(unique(bed$preg.id)),'family_match.svg',sep = "_"), width = 7.5, height = 6)
plot_family_v3(bed)
dev.off()

# 浓度分布图 M-S（V3）
plot_density_MS <- function(prepared_family){
  tmp = copy(prepared_family)[mom.geno=="1/1" & preg.dp >= 10][preg.dp >= quantile(preg.dp, 0.2, na.rm = T) & preg.rf <= quantile(preg.rf, 0.99)]
  tmp = tmp[str_length(mom.geno.bases)==3, .(ref=preg.ref.dp, rf= preg.rf, all=preg.dp)]
  if(length(tmp$rf)>=2) {
  dens = density(tmp$rf)
  dens = data.table(x=dens$x, y=dens$y)
  dens= dens[, d1:= c(0, diff(y))][, d2:= c(diff(y), 0)]
  p = dens[d1 >= 0 & d2 <= 0 & x > 0.005]
  pp = ggplot(tmp) + geom_density(aes(rf))
  if(nrow(p)>=1) {
      pp = pp + geom_text(aes(x=x, y=2*y, label=round(x,4)), data = p) +
        geom_vline(xintercept = p$x, show.legend = T)
    }
  pp
}else{
    pp <- ggplot()
    pp
  }
}
svg(file = paste(sort(unique(bed$dad.id)),sort(unique(bed$mom.id)),sort(unique(bed$preg.id)),'density_MS.svg',sep = "_"), width = 11.78, height = 4)
plot_density_MS(bed)
dev.off()

# 浓度分布图 S-S（V3）
plot_density_SS <- function(prepared_family){
  # 画S-S浓度分布图前，将S的信息赋给M，重新算基因型。M、S基因型计算的差别是S多一条规则，ref、alt都大于err_min时，基因型为0/1,重新计算时去掉这个规则：  preg.bed[ref.dp >= err_min & alt.dp >= err_min,`:=`(geno="0/1", geno.bases=paste0(ref, "/", alt))]
  tmp = copy(prepared_family)
  tmp[,`:=`(mom.ref=preg.ref, mom.alt=preg.alt, mom.dp=preg.dp, mom.ref.dp=preg.ref.dp, mom.alt.dp=preg.alt.dp, mom.rf=preg.rf)]
  tmp[, mom.geno:=""]                                       # 重新计算M基因型前先把M基因型置空(FIX: 线上线下S-S浓度图不一致)，20181113 by hongyuchen
  tmp[mom.dp > 12 & mom.rf >= 0.9, mom.geno:="0/0"]
  tmp[mom.dp >= 12 & mom.rf <= 0.1, mom.geno:="1/1"]
  if(all(c("mom.ref", "mom.alt") %in% names(tmp))){
      tmp[mom.geno=="0/0", mom.geno.bases:=paste0(mom.ref,"/",mom.ref)]
      tmp[mom.geno=="1/1", mom.geno.bases:=paste0(mom.alt,"/",mom.alt)]
     }

  tmp = tmp[mom.geno=="1/1" & preg.dp >= 10][preg.dp >= quantile(preg.dp, 0.2, na.rm = T) & preg.rf <= quantile(preg.rf, 0.99)]
  tmp = tmp[str_length(mom.geno.bases)==3, .(ref=preg.ref.dp, rf= preg.rf, all=preg.dp)]
  if(length(tmp$rf)>=2) {
  dens = density(tmp$rf)
  dens = data.table(x=dens$x, y=dens$y)
  dens= dens[, d1:= c(0, diff(y))][, d2:= c(diff(y), 0)]
  p = dens[d1 >= 0 & d2 <= 0 & x > 0.005]
  pp = ggplot(tmp) + geom_density(aes(rf))
  if(nrow(p)>=1) {
      pp = pp + geom_text(aes(x=x, y=2*y, label=round(x,4)), data = p) +
        geom_vline(xintercept = p$x, show.legend = T)
    }
  pp
}else{
    pp <- ggplot()
    pp
  }
}
svg(file = paste(sort(unique(bed$dad.id)),sort(unique(bed$mom.id)),sort(unique(bed$preg.id)),'density_SS.svg',sep = "_"), width = 11.78, height = 4)
plot_density_SS(bed)
dev.off()

#plot the family
plot_pt <- function(prepared_family, by.dp=F){
  d1 = prepared_family[is.test==T&is.mis==F]
  d2 = prepared_family[is.test==T&is.mis==T]

  p = ggplot(prepared_family[is.test==T], aes(n, abs(dad.rf - preg.rf),
                                              size=mom.dp,
                                              shape=dad.rf>0.95 | dad.rf < 0.05))
  if (nrow(d1)>0) p = p + geom_point(data=d1, colour="green", alpha=I(.5))
  if (nrow(d2)>0) p = p + geom_point(data=d2, colour="red", alpha=I(.5))

  # p + ylim(0, 1) + xlab("") + ylab("") + theme(legend.position="none")
  p + ylim(0, 1)+ xlim(0, 3200) + xlab("") + ylab("") + theme(legend.position="none")  # 修复3200这个标记后面少了个0的bug

}

# 图2
# svg(file=paste(bed$dad.id,bed$mom.id,bed$preg.id,'fig2.svg',sep = "_"))
svg(file=paste(sort(unique(bed$dad.id)),sort(unique(bed$mom.id)),sort(unique(bed$preg.id)),'fig2.svg',sep = "_"))
#cowplot::plot_grid(plot_pt(bed), labels = bed$dad.id[1])
# cowplot::plot_grid(plot_pt(bed), labels = paste(strsplit(bed$dad.id[1], split = '-')[[1]][1:2], collapse = "-"))
cowplot::plot_grid(plot_pt(bed), labels = paste(strsplit(sort(unique(bed$dad.id)), split = '-')[[1]][1:2], collapse = "-"))
dev.off()

# 胎儿浓度
d0 = bed[mom.geno=="1/1" & preg.dp >= 10][preg.dp >= quantile(preg.dp, 0.2, na.rm = T) & preg.rf <= quantile(preg.rf, 0.99)]
d0 = d0[str_length(mom.geno.bases)==3, .(ref=preg.ref.dp, rf= preg.rf, all=preg.dp)]

if(length(d0$rf)>=2) {
dens = density(d0$rf)
dens = data.table(x=dens$x, y=dens$y)
dens= dens[, d1:= c(0, diff(y))][, d2:= c(diff(y), 0)]
p = dens[d1 >= 0 & d2 <= 0 & x > 0.005]

pp = ggplot(d0) + geom_density(aes(rf)) #+geom_point(aes(ref/all, dm*100)) +
if(nrow(p)>=1) {
  pp = pp + geom_text(aes(x=x, y=2*y, label=round(x,4)), data = p) +
    geom_vline(xintercept = p$x, show.legend = T)
}
}else{
  pp <- ggplot()
}
#
# svg(file = paste(bed$dad.id,bed$mom.id,bed$preg.id,'preg_percent.svg',sep = "_"))
svg(file = paste(sort(unique(bed$dad.id)),sort(unique(bed$mom.id)),sort(unique(bed$preg.id)),'preg_percent.svg',sep = "_"))
plot(pp)
dev.off()

info_show <- data.frame(bed$preg.id,dp_preg,percent,error,s_signal,bed$mom.id,dp_mom, mom_preg)

# file_name_txt = paste(bed$mom.id,bed$preg.id,'info_show.txt',sep = "_")
file_name_txt = paste(sort(unique(bed$mom.id)),sort(unique(bed$preg.id)),'info_show.txt',sep = "_")
write.table(info_show[c(2),], file = file_name_txt, row.names = F, quote = F, sep="\t")


trim_geno <- function(geno) {
  stringr::str_replace_all(geno, ",.[^/]*", "")
}

format_family <- function(bed, rs.included=F) {
  out = bed[order(n),][is.test==T, .(`检测位点编号`=n,
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
file_name = paste(sort(unique(bed$dad.id)),sort(unique(bed$mom.id)),sort(unique(bed$preg.id)),'test_pos.txt',sep = "_")
write.table(test_pos, file = file_name, row.names = F, quote = F, sep="\t")
# file_name = paste(bed$dad.id,bed$mom.id,bed$preg.id,'test_pos.txt',sep = "_")
# write.table(test_pos, file = file_name[1], row.names = F, quote = F, sep="\t")
