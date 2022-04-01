library(RTCGA.rnaseq)
library(dplyr)

# 构成20533个基因在三类癌症中的测序表达谱
extract.list <- expressionsTCGA(BRCA.rnaseq, OV.rnaseq, HNSC.rnaseq)

# 查看三类癌症中的样本数量
n.samples <- table(extract.list$dataset)

# 将extract.list对象中的"dataset"向量重命名为"cohort"向量
extract.list <- dplyr::rename(extract.list, cohort = dataset)

# 找到extract.list对象中的"bcr_patient_barcode"向量中坐标于14-15的字符串等于"01"的样本（编号01~09表示肿瘤，10~19表示正常对照）
cancer.list <- filter(extract.list, substr(bcr_patient_barcode, 14, 15) == "01")

# 进行主成分分析，使用不同的颜色标记不同的癌症
pca.plot <- pcaTCGA(cancer.list, "cohort")
plot(pca_plot)
