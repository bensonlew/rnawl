# use: pt_dup.R WQ123M.tab WQ123S.tab 2 /mnt/ilustre/users/sanger-dev/sg-users/zhoumoli/pt/targets.bed.rda
# /mnt/ilustre/users/sanger-dev/sg-users/zhoumoli/pt/targets.bed.rda

library(RcppRoll)
library(magrittr)
library(stringr)
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(dtplyr)

args <- commandArgs(TRUE)
mom.tab <- args[1]
cat(mom.tab)
cat("\t")
preg.tab <- args[2]
cat(preg.tab)
cat("\t")
err.min <- as.numeric(args[3])
cat(err.min)
cat("\t")
load(args[4])
cat(args[4])
cat("\t")
out_path = args[5]
dir.create(out_path)
cat(out_path)
file_path = args[6]
cat(file_path)
cat("\n")
dad_id = args[7]
cat(dad_id)


#file_list <- list.files('/mnt/ilustre/users/sanger-dev/app/database/human/pt_ref/tab_data/') # 测试机
#file_list <- list.files('/mnt/ilustre/users/sanger/app/database/human/pt_ref/tab_data/') # 正式机
file_list <- list.files(args[6])
dad.ids <- character(0)
for (m in file_list){dad.ids <- append(dad.ids, strsplit(m, split="\\.")[[1]][1])} 
get_bed = function(.id, dir=args[6]){ paste0(dir, .id, ".tab")}
#get_bed = function(.id, dir="/mnt/ilustre/users/sanger-dev/app/database/human/pt_ref/tab_data/"){ paste0(dir, .id, ".tab")} # 测试机
#get_bed = function(.id, dir="/mnt/ilustre/users/sanger/app/database/human/pt_ref/tab_data/"){ paste0(dir, .id, ".tab")} # 正式机
#' read bedfile
#' @export
read_bed <- function(.bedFile){
  .dt = fread(.bedFile, sep="\t") %>%
    setnames(paste0("V",1:8),c("id","chrom","pos","ref","alt","dp","ref.dp","alt.dp"))

  .dt$chrom = as.character(.dt$chrom)
  .dt
}

add_geno <- function(.bed, err.min=2, rf.min=0.1, rf.max=0.9, dp.min=12){
  if(nrow(.bed)<=0)  {
    stop("nrow <=0")
  }
  if(!all(c("chrom","pos","ref.dp","alt.dp") %in% names(.bed))){
    stop("must have chrom,post,ref.dp,alt.dp columns.")
  }
  .bed[alt==".", alt:=ref]
  .bed[, dp:=alt.dp+ref.dp]
  .bed[, rf:=ref.dp/dp]
  .bed[, geno:=paste0(ifelse(ref.dp >= err.min, 0, "*"),
                      "/",
                      ifelse(alt.dp >= err.min, 1, "*"))]
  .bed[dp > dp.min & rf >= rf.max, geno:="0/0"]
  .bed[dp >= dp.min & rf <= rf.min, geno:="1/1"]
  if(all(c("ref", "alt") %in% names(.bed))){
    .bed[dp <= dp.min, geno.bases:=paste0(ifelse(ref.dp >=err.min, ref,"*"),
                                          "/",
                                          ifelse(alt.dp>=err.min, alt, "*"))]
    .bed[geno=="0/0", geno.bases:=paste0(ref,"/",ref)]
    .bed[geno=="0/1", geno.bases:=paste0(ref,"/",alt)]
    .bed[geno=="1/1", geno.bases:=paste0(alt,"/",alt)]
  }
  .bed
}

annotate_family <- function(.bed) {
  .bed = .bed[targets.bed, on=c("chrom","pos")]
  .bed[, `:=`(is.test=F, is.mis=F, mustbe=F, mustnotbe=F, good=F)]
  .bed[mom.geno %in% c("0/0", "1/1") &
         preg.geno %in% c("0/1") &
         dad.dp > 15 &
         (dad.geno %in% c("0/0") | (dad.geno %in% c("0/1", "1/1") & str_sub(dad.alt,1,1) == str_sub(preg.alt,1,1))) &
         mj.rf >=0.1 & mj.rf <= .9 &
         (preg.geno.bases == hapmap.geno | is.na(hapmap.geno)) &
         !str_detect(dad.geno, "[*]") &
         !str_detect(preg.geno.bases, ",") &
         !chrom %in% c("chrX", "chrY") &
         str_length(trim_geno(preg.ref))==1 & str_length(trim_geno(preg.alt))==1 &
         str_length(trim_geno(mom.ref))==1 & str_length(trim_geno(mom.alt))==1 &
         str_length(trim_geno(dad.ref))==1 & str_length(trim_geno(dad.alt))==1
       , is.test:= T]


  .bed[is.test & mom.geno=="0/0" & dad.geno == "0/0", `:=`(pi=0.05, is.mis=T)]
  .bed[is.test & mom.geno=="0/0" & dad.geno == "0/1", `:=`(pi=0.5/(1-mj.rf))]
  .bed[is.test & mom.geno=="0/0" & dad.geno == "1/1", `:=`(pi=1/(1-mj.rf), good=T)]
  .bed[is.test & mom.geno=="1/1" & dad.geno == "0/0", `:=`(pi=1/mj.rf, good=T)]
  .bed[is.test & mom.geno=="1/1" & dad.geno == "0/1", `:=`(pi=0.5/mj.rf)]
  .bed[is.test & mom.geno=="1/1" & dad.geno == "1/1", `:=`(pi=0.05, is.mis=T)]

  .bed[preg.dp > 30 & mom.geno=="0/0" & dad.geno == "1/1", `:=`(mustbe=T)]
  .bed[preg.dp > 30 & mom.geno=="1/1" & dad.geno == "0/0", `:=`(mustbe=T)]

  .bed[preg.dp > 30 & mom.geno=="0/0" & dad.geno == "0/0", `:=`(mustnotbe=T)]
  .bed[preg.dp > 30 & mom.geno=="1/1" & dad.geno == "1/1", `:=`(mustnotbe=T)]

  .bed
}

#' trim genotype
#' @export
trim_geno <- function(geno) {
  stringr::str_replace_all(geno, ",.[^/]*", "")
}

join_family <- function(dad.bed, mom.bed, preg.bed, err.min=2){
  dad.bed = copy(dad.bed)
  mom.bed = copy(mom.bed)
  preg.bed = copy(preg.bed)

  #correct_geno(preg.bed)
  preg.bed[ref.dp >= err.min & alt.dp >= err.min,
           `:=`(geno="0/1", geno.bases=paste0(ref, "/", alt))]

  .family = list(dad=dad.bed, mom=mom.bed, preg=preg.bed)

  # set names
  by.names = c("chrom", "pos")
  for(mem in c("dad", "mom", "preg")) {
    mem.names = setdiff(names(.family[[mem]]), by.names)
    setnames(.family[[mem]], mem.names, paste0(mem, ".", mem.names))
  }

  joined.family = full_join(.family$dad, .family$preg, by=by.names) %>%
    full_join(.family$mom, by=by.names)

  as.data.table(joined.family)[, reg:=paste0(chrom,":",pos)]
}

# for (i in 1:length(dad.ids)){
  # filename = paste0("/mnt/ilustre/users/sanger-dev/app/database/human/pt_ref/tab_data/",dad.ids[i],".tab")
  # if (file.info(filename)$size == FALSE ){
    # print(dad.ids[i])
    # print("the tab is empty")
    # dad.ids <- dad.ids[-i]
  # }
# }

mom.bed =  mom.tab %>% read_bed() %>% add_geno()
preg.bed = preg.tab %>% read_bed() %>% add_geno(err.min=err.min)

dup.family = dad.ids %>% lapply(function(dad.id){
  dad.bed = dad.id %>% get_bed() %>% read_bed() %>% add_geno()
  join_family(dad.bed, mom.bed, preg.bed, err.min=err.min) %>% annotate_family()
}) %>% rbindlist(fill=T) %>% as.data.table()

output <- data.frame(dup.family[!is.na(dad.id), .(
  # dad.name = samples[gid==dad.id,name],
  # dad.name = dad.id
  test.pos.n=sum(is.test),
  err.pos.n=sum(is.mis),
  cpi=prod(pi, na.rm = T),
  dp=mean(dad.dp, na.rm = T),
  good=paste0(round(sum(is.test & mustbe,na.rm = T)/sum(mustbe,na.rm = T)*100, 2), "% : ",sum(is.test & mustbe, na.rm = T),"/",sum(mustbe, na.rm = T)),
  bad=paste0(round(sum(is.test&mustnotbe,na.rm = T)/sum(mustnotbe,na.rm = T)*100, 2), "% : ",sum(is.test & mustnotbe, na.rm = T),"/",sum(mustnotbe, na.rm = T))),
  by=.(dad.id)
  ][, .(
    #`父本名称`="dad.name",
    `父本ID`=dad.id,
    `测试位点数`=test.pos.n,
    `错配位点数`=err.pos.n,
    `错配率(%)`=(err.pos.n / test.pos.n*100)%>% round(2),
    `父权值`=(cpi),
    `深度`=dp%>%round(2),
    `有效率`=good,
    `无效率`=bad,
    `匹配`=ifelse(err.pos.n / test.pos.n < 0.15 & (cpi) > 10000, "Yes","No"))
    ][order(`错配率(%)`, -`父权值`)]
)

mom_name <- strsplit(strsplit(mom.tab, split='/')[[1]][length(strsplit(mom.tab, split='/')[[1]])], split='\\.')[[1]][1]
preg_name <- strsplit(strsplit(preg.tab, split='/')[[1]][length(strsplit(preg.tab, split='/')[[1]])], split='\\.')[[1]][1]
outputfile_name = paste0("./", out_path, "/", dad_id, "_", mom_name, "_", preg_name, ".txt")
write.table(output,file = (outputfile_name),sep = "\t", row.names = F, quote = F)
