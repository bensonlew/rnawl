library(magrittr)
library(data.table,warn.conflicts = FALSE)
library(RcppRoll) #移动平均��?
library(stringr)  #处理字符��?
library(dtplyr)
library(dplyr,warn.conflicts = FALSE)    #处理表格
library(ggplot2)
library(cowplot,warn.conflicts = FALSE)
library(Cairo)
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
  d.dads = joined_family[dad.geno %in% c("0/0","1/1", "0/1")][, .(dp=dad.dp, n, rf=dad.rf, name=paste(strsplit(dad.id, split = '-')[[1]][1:2], collapse = "-"), alpha=1)]
  #View(d.dads)
  d.mom = joined_family[mom.geno %in% c("0/0","1/1", "0/1"), .(dp= mom.dp, n, rf=mom.rf, name=paste(strsplit(mom.id, split = '-')[[1]][1:2], collapse = "-"), alpha=1)]
  #View(d.mom)
  d.preg = joined_family[, .(dp= preg.dp, n, rf=preg.rf, name=paste(strsplit(preg.id, split = '-')[[1]][1:2], collapse = "-"), alpha=1+is.test)]
  #View(d.preg)
  if (by.dp)
    p <- qplot(dp, rf, data=rbind(d.dads, d.mom, d.preg)[!is.na(name)], facets = name ~ ., alpha=I(.1)) + geom_point(aes(size= alpha, colour=alpha), alpha=I(.2))
  else
    p <- qplot(n, rf, data=rbind(d.dads, d.mom, d.preg)[!is.na(name)], facets = name ~ ., alpha=I(.1)) + geom_point(aes(size= alpha, colour=alpha), alpha=I(.2))
  p + xlab("") + ylab("") + theme(legend.position="none")
}

# 家系��?
svg(file = paste(bed$dad.id,bed$mom.id,bed$preg.id,'family.svg',sep = "_"))
plot_family(bed, by.dp = T)
dev.off()

#图1
svg(file=paste(bed$dad.id,bed$mom.id,bed$preg.id,'fig1.svg',sep = "_"))
plot_family(bed, by.dp = F)
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

  p + ylim(0, 1) + xlab("") + ylab("") + theme(legend.position="none")
}

# 图2
svg(file=paste(bed$dad.id,bed$mom.id,bed$preg.id,'fig2.svg',sep = "_"))
#cowplot::plot_grid(plot_pt(bed), labels = bed$dad.id[1])
cowplot::plot_grid(plot_pt(bed), labels = paste(strsplit(bed$dad.id[1], split = '-')[[1]][1:2], collapse = "-"))
dev.off()

# 胎儿浓度
d0 = bed[mom.geno=="1/1" & preg.dp >= 10][preg.dp >= quantile(preg.dp, 0.2, na.rm = T) & preg.rf <= quantile(preg.rf, 0.99)]
d0 = d0[str_length(mom.geno.bases)==3, .(ref=preg.ref.dp, rf= preg.rf, all=preg.dp)]

dens = density(d0$rf)
dens = data.table(x=dens$x, y=dens$y)
dens= dens[, d1:= c(0, diff(y))][, d2:= c(diff(y), 0)]
p = dens[d1 >= 0 & d2 <= 0 & x > 0.005]

pp = ggplot(d0) + geom_density(aes(rf)) #+geom_point(aes(ref/all, dm*100)) +
if(nrow(p)>=1) {
  pp = pp + geom_text(aes(x=x, y=2*y, label=round(x,4)), data = p) +
    geom_vline(xintercept = p$x, show.legend = T)
}
#
svg(file = paste(bed$dad.id,bed$mom.id,bed$preg.id,'preg_percent.svg',sep = "_"))
plot(pp)
dev.off()

info_show <- data.frame(bed$preg.id,dp_preg,percent,error,s_signal,bed$mom.id,dp_mom, mom_preg)

file_name_txt = paste(bed$mom.id,bed$preg.id,'info_show.txt',sep = "_")
write.table(info_show[c(2),], file = file_name_txt[1], row.names = F, quote = F, sep="\t")


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
file_name = paste(bed$dad.id,bed$mom.id,bed$preg.id,'test_pos.txt',sep = "_")
write.table(test_pos, file = file_name[1], row.names = F, quote = F, sep="\t")
