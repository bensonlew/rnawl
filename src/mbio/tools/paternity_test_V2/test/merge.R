options(stringsAsFactors=F)

#args <- commandArgs(T) #获取参数，此时args[1]是up,args[2]是down
up_pd <- read.csv("20170924.p.csv",header=T,check.names=F,fileEncoding="GBK")
down_pd <- read.csv("20170924.t.csv",header=T,check.names=F,fileEncoding="GBK",skip=2)
#up_pd <- read.csv(args[1],header=T,check.names=F,fileEncoding="GBK")
#down_pd <- read.csv(args[2],header=F,check.names=F,fileEncoding="GBK",skip=2)

##删除NA行
down_pd <- down_pd[apply(down_pd,1,function(x)x[1] != ""),]
down_pd <- down_pd[,1:12]
down_pd <- as.matrix(cbind(down_pd[,c(6,6,7)],down_pd))
up_pd <- as.matrix(up_pd[,1:15])
interval_file <- rbind(up_pd,down_pd)
colnames(interval_file)[(dim(interval_file)[2]-1):dim(interval_file)[2]] = c("分析类型","备注")
interval_file <- as.data.frame(interval_file)
rownames(interval_file) <- 1:dim(interval_file)[1]

###开始做数据检查
##1. 检查有没有相同的index
table(interval_file$'index序列')