#for example: Rscript "1,2,3,4,5" 1 /mnt/ilustre/users/sanger-dev/app/database/human/hg38_nipt/bed_file /mnt/ilustre/users/sanger-dev/app/database/human/hg38_nipt/bed_file
library(RcppRoll)  ##利用滑动窗口（Efficient windowed / rolling operations）
library(magrittr)  ##一个高效的管道操作工具包
library(stringr)  ##一个字符串处理工具集
library(dplyr)  ##提供了一个访问常见数据库的接口
library(data.table) 

args <- commandArgs(T)
test.gids <- args[1]  # 样本名字
ref_group <- args[2]  # 第几个参考组，以及用于设定参考组的名字
bed.dir <- args[3] # 样本存放的地址
ref.dir <- args[4] # 参考组所在目录
cat(test.gids)
cat("\n")
cat(ref_group)
cat("\n")
cat(bed.dir)
cat("\n")
cat(ref.dir)
cat("\n")

get.bed = function(gid) file.path(bed.dir, paste0(gid,".bed.2"))  ##获取样本的bed文件函数

ref.cor.cn.file = file.path(ref.dir, paste0(ref_group, ".ref.cor.Rdata"))

##把bed文件的每一列分割处理
fread.bed = function(file) {
  message("fread ", file)
  x = fread(file, header = F) %>%
    setnames(c("V1","V2","V3","V4","V5","V6","V7","V8"),
             c("chr","start","end","gc","map","pn","reads","gid"))
  x$gid = as.character(x$gid)
  x
}
##计算百分比
cor.cn = function (x, map.min = 0.9, samplesize = 50000,
                   verbose = T, gc.cor = T, map.cor = T) {
  message("run cor.cn for ", unique(x$gid))

  if (length(x$reads) == 0 | length(x$gc) == 0 | length(x$map) == 0)
    stop("Missing one of required columns: reads, gc, map")
  x = data.table(x)
  x[, valid:=(reads>0 & gc > 0)]

  reads.range <- x[, quantile(reads[valid], prob = c(0.01, 1-0.01), na.rm = TRUE)]
  gc.range <- x[, quantile(gc[valid], prob = c(0.01, 1 - 0.01), na.rm = TRUE)]
  x[, ideal:=(map > map.min &
               reads > reads.range[1] & reads < reads.range[2] &
               chr <= 22 &
               gc > gc.range[1] & gc < gc.range[2] &
               pn == 0)]

  set <- which(x$ideal)
  select <- sample(set, min(length(set), samplesize))
  rough <- loess(x$reads[select] ~ x$gc[select], span = 0.03)
  i <- seq(0, 1, by = 0.001)
  final <- loess(predict(rough, i) ~ i, span = 0.3)
  x[reads != 0 & gc != 0, cor.gc := reads/predict(final, gc)]

  coutlier <- 0.001
  range <- quantile(x$cor.gc[which(x$ideal)], prob = c(coutlier, 1 - coutlier), na.rm = TRUE)
  set <- which(x$cor.gc < range[2] & x$cor.gc > range[1])
  select <- sample(set, min(length(set), samplesize))
  final <- approxfun(lowess(x$map[select], x$cor.gc[select]))
  x[, cor.map := cor.gc/final(map)][,copy:=map.cor*cor.map+(!map.cor)*cor.gc]
}
##进行划窗
bin.cn <- function(x, bin.step=1*10^6, bin.window=bin.step){
  message("run bin.cn")
  x = data.table(x)[map>0.8,]
  x = x[, bin:= start %/% bin.step+1][
    , .(cn1=mean(copy, na.rm = T)), by=.(chr,bin)]
  cn1.mean = x[chr %in% c(1:12,14:17,19:20,22),][cn1 > quantile(cn1,0.05,na.rm = T) & cn1 < quantile(cn1,0.95,na.rm = T),mean(cn1)]
  x[,cn1:=cn1/cn1.mean]

  n <- bin.window %/% bin.step
  x.1 = x[chr!=25][
    , .(cn=roll_mean(cn1, n, na.rm = T)), by=.(chr)][
    , bin:=seq_along(cn), by=.(chr)][
    , n:=seq_along(cn),
   ]
  x.2 = x[chr==25][
      , .(cn=cn1), by=.(chr)][
      , bin:=seq_along(cn), by=.(chr)][
      , n:=seq_along(cn)+max(x.1$n),
   ]
  list(x.1, x.2) %>% rbindlist(use.names=T)
}

ref.gids <- strsplit(test.gids, split=',')[[1]]
message("start update ref")
system.time({x = ref.gids %>% lapply(get.bed) %>% lapply(fread.bed) %>% lapply(cor.cn) %>% rbindlist()})
# x = ref.gids %>% lapply(function(ref.id){x = ref.id %>% get.bed() %>% fread.bed() %>% cor.cn()}) %>% rbindlist()
message("end")
save(x, file=ref.cor.cn.file)
message("save success")




