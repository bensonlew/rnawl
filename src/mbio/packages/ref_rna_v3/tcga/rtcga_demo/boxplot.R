library(RTCGA)
library(RTCGA.mRNA)
library(RTCGA.rnaseq)
library(ggpubr)

# 构成五个基因在三类癌症中的芯片表达谱
expr <- expressionsTCGA(BRCA.mRNA, OV.mRNA, LUSC.mRNA,
                        extract.cols = c("GATA3", "PTEN", "XBP1", "ESR1", "MUC1"))

# 查看三类癌症中的样本数量
nb_samples <- table(expr$dataset)

# 对"dataset"和"bcr_patient_barcode"两列重命名
expr$dataset <- gsub(pattern = ".mRNA", replacement = "", expr$dataset)
expr$bcr_patient_barcode <- paste0(expr$dataset, c(1:590, 1:561, 1:154))

# 箱型图："GATA3"在三类癌症中的表达量分布
ggboxplot(expr, x = "dataset", y = "GATA3", title = "GATA3", ylab = "Expression", color = "dataset", palette = "jco")

# 箱型图："PTEN"在三类癌症中的表达量分布
ggboxplot(expr, x = "dataset", y = "PTEN", title = "PTEN", ylab = "Expression", color = "dataset", palette = "npg")

# 箱型图："XBP1"在三类癌症中的表达量分布
ggboxplot(expr, x = "dataset", y = "XBP1", title = "XBP1", ylab = "Expression", color = "dataset", palette = "aaas")

# 箱型图："ESR1"在三类癌症中的表达量分布
ggboxplot(expr, x = "dataset", y = "ESR1", title = "ESR1", ylab = "Expression", color = "dataset", palette = "lancet")

# 箱型图："MUC1"在三类癌症中的表达量分布，展示比较的P值
my_comparisons <- list(c("BRCA", "OV"), c("OV", "LUSC"))
ggboxplot(expr, x = "dataset", y = "MUC1", title = "MUC1", ylab = "Expression", color = "dataset", palette = "ucscgb") + stat_compare_means(comparisons = my_comparisons)

# 选定基因，进行不同癌症之间的均值比较
df_cm <- compare_means(c(GATA3, PTEN, XBP1) ~ dataset, data = expr)

# 分面箱型图："GATA3"、"PTEN"和"XBP1"三个基因的表达量分布，标记出在"BRCA"和"OV"中表达量大于3.9的样本
label.select.criteria <- list(criteria = "`y` > 3.9 & `x` %in% c('BRCA', 'OV')")
ggboxplot(expr, x = "dataset",
          y = c("GATA3", "PTEN", "XBP1"),
          combine = TRUE,
          color = "dataset", palette = "jco",
          ylab = "Expression",
          label = "bcr_patient_barcode", # column containing point labels
          label.select = label.select.criteria, # Select some labels to display
          font.label = list(size = 9, face = "italic"), # label font
          repel = TRUE # Avoid label text overplotting
)

# 构成五个基因在三类癌症中的测序表达谱
seq_expr <- expressionsTCGA(BRCA.rnaseq, OV.rnaseq, LUSC.rnaseq,
                            extract.cols = c("GATA3|2625", "PTEN|5728", "XBP1|7494","ESR1|2099", "MUC1|4582"))

# 查看三类癌症中的样本数量
nb_seq_samples <- table(seq_expr$dataset)

# 箱型图："PTEN|5728"在三类癌症中的表达量分布
ggboxplot(seq_expr, x = "dataset", y = "`PTEN|5728`", title = "ESR1|2099", ylab = "Expression", color = "dataset", palette = "jco")
