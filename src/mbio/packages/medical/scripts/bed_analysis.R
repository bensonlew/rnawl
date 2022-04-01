################该脚本用于对nipt的bed文件进行统计分析########################
#example：Rscript bed_analysis.R D:\\majorbio\\nipt_test_data\\150_R1.bed.2 D:\\majorbio\\nipt_test_data\\2.ref.cor.Rdata 10 1 nipt_output

library(RcppRoll)  ##利用滑动窗口（Efficient windowed / rolling operations）
library(magrittr)  ##一个高效的管道操作工具包
library(stringr)  ##一个字符串处理工具集
library(dplyr)  ##提供了一个访问常见数据库的接口
library(data.table)  ##用于快速处理大数据集;fread函数读取数据快，处理数据框也快
args <- commandArgs(T)
test.gids <- args[1]  # 样本的绝对路径
ref.cor.cn.file <- args[2]  #参考组的文件路径
bw <- as.numeric(args[3])   #滑窗大小
bs <- as.numeric(args[4])   #滑窗步长
outdir <- args[5] # 输出目录
dir.create(outdir)
cat(args[1])
cat("\n")
cat(args[2])
cat("\n")
cat(args[3])
cat("\n")
cat(args[4])
cat("\n")
cat(args[5])
load(ref.cor.cn.file) 
# get.bed = function(gid) file.path(bed.dir, paste0(gid,".bed.2"))  ##获取样本的bed文件函数
 
#这一步是用于更新参考组的 
# ref.cor.cn = reactive({
# ref.cor.cn.file = file.path(bed.dir, "..", "data", paste0(input$ref_group, ".ref.cor.Rdata"))
#     load(ref.cor.cn.file) 
#     ##这个判断应该是单独写成一个tool
#     if (!setequal(as.character(ref.gids()),x$gid%>% unique %>% as.character)) {
#       message("updating ref")
#       x = ref.gids() %>%
#         lapply(get.bed) %>% lapply(fread.bed) %>%
#         lapply(cor.cn) %>%
#         rbindlist()
#       save(x, file=ref.cor.cn.file)
#     }
#     x
#   })


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
##计z-score值
z.cn <- function(test.cn, ref.cn, bs, bw){
  message("run z.cn")
  message("ref.gids = ", paste(ref.cn$gid %>% unique(), collapse=","))
  message("test.gids = ", paste(test.cn$gid %>% unique(), collapse=","))
  ref.binned.cns = data.table(ref.cn)[, bin.cn(.SD,bs,bw), by=gid]
  test.binned.cns = data.table(test.cn)[, bin.cn(.SD,bs,bw), by=gid]

  ref = ref.binned.cns[, .(sd=sd(cn,na.rm = T), mean=mean(cn, na.rm = T)),
                       by=.(chr, bin)]
  z = test.binned.cns %>%
    left_join(ref, by=c("chr","bin")) %>%
    group_by(gid) %>%
    mutate(z=(cn-mean)/sd,n=seq_along(cn))
  z
}

test.cor.cn = test.gids %>% lapply(fread.bed) %>% lapply(cor.cn) %>% rbindlist()
##这个就是画图结果表zs.out
zs.out <- z.cn(test.cor.cn, x, bs=bs*10^6, bw=bw*10^6)
write.table(zs.out, file = paste(outdir, "/z.xls", sep = ""), row.names = F, quote = F, sep = "\t")
##分析结果（表）,就是zs.out
##对应的是网站中统计
zs = zs.out %>% as.data.table()
zz_result <- zs[chr <= 22, .(zz = sqrt(mean(z^2))), by=.(gid)]
write.table(zz_result, file = paste(outdir, "/zz.xls", sep = ""), row.names = F, quote = F, sep = "\t")

##拷贝数变异（图）, js画