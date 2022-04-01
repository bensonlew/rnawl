library(RTCGAToolbox)

# 载入下载完毕的"BRCA"队列的RNA测序数据、临床数据、测序突变数据和处理过的拷贝数数据
BRCA.data <- getFirehoseData("BRCA", RNASeqGene = TRUE, clinical = TRUE, Mutation = TRUE, GISTIC = TRUE)

# 提取测序突变数据（）
Mutation.data <- biocExtract(BRCA.data, "Mutation")

getCNGECorrelation(BRCA.data)

# 查看给定癌症队列中每一个基因的突变率
mutRate <- getMutationRate(BRCA.data)


