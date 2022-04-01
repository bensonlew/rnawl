library(magrittr)
library(data.table,warn.conflicts = FALSE)
library(RcppRoll) #移动平均��?
library(stringr)  #处理字符��?
library(dtplyr)
library(dplyr,warn.conflicts = FALSE)    #处理表格
library(ggplot2)
library(cowplot,warn.conflicts = FALSE)
library(xlsx)
# /mnt/ilustre/users/sanger-dev/app/program/R-3.3.1/bin/Rscript family_joined.R /mnt/ilustre/users/sanger-dev/workspace/20161219/PtProcess_pt_process_test/output/WQ2130-FC.tab /mnt/ilustre/users/sanger-dev/workspace/20161219/PtProcess_pt_process_test/output/WQ2130-M.tab /mnt/ilustre/users/sanger-dev/workspace/20161219/PtProcess_pt_process_test/output/WQ2130-S.tab 2 /mnt/ilustre/users/sanger-dev/sg-users/zhoumoli/pt/targets.bed.rda
# modified by hongdong 20180106 处理了父本母本深度低的的时候，基因分型不正确的bug
args <- commandArgs()
file_F = args[6]
cat("\n")
file_M = args[7]
cat("\n")
file_S = args[8]
cat("\n")
err_min = as.numeric(args[9])
rda_file = args[10]
#out_path = args[11]

load(args[10])
trim_geno <- function(geno) {
  stringr::str_replace_all(geno, ",.[^/]*", "")
}

# 更改表头
get_bed <- function(file, err_min=err_min, rf.min=0.1, rf.max=0.9, dp.min=12){
    dt = fread(file,sep="\t") %>% setnames(paste0("V",1:8),c("id","chrom","pos","ref","alt","dp","ref.dp","alt.dp"))
    dt$chrom = as.character(dt$chrom)

    if(nrow(dt)<=0){
      stop("nrow <= 0")
    }
    if(!all(c("chrom","pos","ref.dp","alt.dp") %in% names(dt))){
      stop("must have chrom,post,ref.dp,alt.dp columns.")
    }

    dt[alt==".",alt := ref] #如果测序测出来的位点是空，就用ref补起��?
    dt[,dp:=alt.dp+ref.dp]# 更改总测序深度为ref与alt测序深度相加
    dt[,rf:=ref.dp/dp]
    dt[,geno:=paste0(ifelse(ref.dp >= err_min,0,"*"),
                     "/",
                     ifelse(alt.dp >= err_min,1,"*"))]
    dt[dp > 12 & rf >= 0.9, geno:="0/0"]
    dt[dp >= 12 & rf <= 0.1, geno:="1/1"]
    if(all(c("ref", "alt") %in% names(dt))){
      dt[dp <= 12, geno.bases:=paste0(ifelse(ref.dp >=err_min, ref,"*"),
                                            "/",
                                            ifelse(alt.dp>=err_min, alt, "*"))]
      dt[geno=="0/0", geno.bases:=paste0(ref,"/",ref)]
      dt[geno=="0/1", geno.bases:=paste0(ref,"/",alt)]
      dt[geno=="1/1", geno.bases:=paste0(alt,"/",alt)]
     }
}

#######joined family函数
dad.bed <- get_bed(file_F, err_min=2)
mom.bed <- get_bed(file_M, err_min=2)
preg.bed <- get_bed(file_S, err_min=err_min)

dad.bed = copy(dad.bed)
mom.bed = copy(mom.bed)
preg.bed = copy(preg.bed)

preg.bed[ref.dp >= err_min & alt.dp >= err_min,`:=`(geno="0/1", geno.bases=paste0(ref, "/", alt))]

family = list(dad=dad.bed, mom=mom.bed, preg=preg.bed)
by.names = c("chrom", "pos")
for(mem in c("dad", "mom", "preg")) {
    mem.names = setdiff(names(family[[mem]]), by.names)
    setnames(family[[mem]], mem.names, paste0(mem, ".", mem.names))
}
joined.family = full_join(family$dad, family$preg, by=by.names) %>%
  full_join(family$mom, by=by.names)
  as.data.table(joined.family)
joined.family[, reg:=paste0(chrom,":",pos)] #上线要改

########annotate family函数
bed = as.data.table(joined.family)
bed = bed[targets.bed, on=c("chrom","pos")]
bed[,`:=`(is.test=F, is.mis=F,mustbe=F,mustnotbe=F, good=F)]

###符合要求的位点标记test为True
bed[mom.geno %in% c("0/0","1/1") & #母本纯合
   preg.geno %in% c("0/1") & #胎儿杂合
   dad.dp > 15 &
   abs(mom.rf-preg.rf) <= 0.25 &
   (dad.geno %in% c("0/0") | (dad.geno %in% c("0/1","1/1") &
   str_sub(dad.alt,1,1) == str_sub(preg.alt,1,1))) &
   mj.rf >=0.1 & mj.rf <= .9 &
   (preg.geno.bases == hapmap.geno | is.na(hapmap.geno)) &
   !str_detect(dad.geno, "[*]") &
   !str_detect(preg.geno.bases, ",") &
   !chrom %in% c("chrX", "chrY") &
   str_length(trim_geno(preg.ref))==1 & str_length(trim_geno(preg.alt))==1 &
   str_length(trim_geno(mom.ref))==1 & str_length(trim_geno(mom.alt))==1 &
   str_length(trim_geno(dad.ref))==1 & str_length(trim_geno(dad.alt))==1
 , is.test:= T
]

bed[is.test & mom.geno=="0/0" & dad.geno=="0/0", `:=`(pi=0.05, is.mis = T)]
bed[is.test & mom.geno=="0/0" & dad.geno=="0/1",`:=`(pi=0.5/(1-mj.rf))]
bed[is.test & mom.geno=="0/0" & dad.geno=="1/1",`:=`(pi=1/(1-mj.rf),good=T)]
bed[is.test & mom.geno=="1/1" & dad.geno == "0/0", `:=`(pi=1/mj.rf, good=T)]
bed[is.test & mom.geno=="1/1" & dad.geno == "0/1", `:=`(pi=0.5/mj.rf)]
bed[is.test & mom.geno=="1/1" & dad.geno == "1/1", `:=`(pi=0.05, is.mis=T)]
#mustbe 是为了有效率，mustbotbe 是为了计算无效率
bed[preg.dp > 30 & mom.geno=="0/0" & dad.geno == "1/1", `:=`(mustbe=T)]
bed[preg.dp > 30 & mom.geno=="1/1" & dad.geno == "0/0", `:=`(mustbe=T)]

bed[preg.dp > 30 & mom.geno=="0/0" & dad.geno == "0/0", `:=`(mustnotbe=T)]
bed[preg.dp > 30 & mom.geno=="1/1" & dad.geno == "1/1", `:=`(mustnotbe=T)]


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
info_show <- data.frame(bed$preg.id,dp_preg,percent,error,s_signal,bed$mom.id,dp_mom, mom_preg)
file_name_txt = paste(sort(unique(bed$mom.id)),sort(unique(bed$preg.id)),'info_show.txt',sep = "_")
write.table(info_show[c(2),], file = file_name_txt, row.names = F, quote = F, sep="\t")
