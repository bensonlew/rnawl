source("https://bioconductor.org/biocLite.R")
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

biocLite("RTCGAToolbox")

library(RTCGAToolbox)

# 函数返回一个包含癌症队列的字符向量
firehose.datasets <- getFirehoseDatasets()

# 函数返回标准数据发布日期的字符向量
firehose.running.dates <- getFirehoseRunningDates()

# 函数返回分析发布日期的字符向量
firehose.analyze.dates <- getFirehoseAnalyzeDates()

# 下载"BRCA"队列的临床数据和测序突变数据
BRCA.data <- getFirehoseData("BRCA", RNASeqGene = TRUE, clinical = TRUE, Mutation = TRUE, GISTIC = TRUE)

# 提取临床数据
clinical.data <- biocExtract(BRCA.data, "clinical")

# 提取测序突变数据
Mutation.data <- biocExtract(BRCA.data, "Mutation")
