library(dplyr)
library(RTCGA.rnaseq)

# 整理3个基因在4种癌症中的表达谱，并将一个"MET|4233"表达量以区间展示
expressionsTCGA(ACC.rnaseq, BLCA.rnaseq, BRCA.rnaseq, OV.rnaseq,
                extract.cols = c("MET|4233", "ZNF500|26048", "ZNF501|115560")) %>%
    rename(cohort = dataset) %>%
    filter(substr(bcr_patient_barcode, 14, 15) == "01") %>%
    mutate(MET = cut(`MET|4233`, round(quantile(`MET|4233`, probs = seq(0, 1, 0.25)), -2),
                     include.lowest = TRUE, dig.lab = 5)) -> ACC_BLCA_BRCA_OV.rnaseq

# 按癌种类型和"MET|4233"表达区间分组，在每个组中对"ZNF500|26048"和"ZNF501|115560"求均值
ACC_BLCA_BRCA_OV.rnaseq %>%
    select(-bcr_patient_barcode) %>%
    group_by(cohort, MET) %>%
    summarise_each(funs(median)) %>%
    mutate(ZNF500 = round(`ZNF500|26048`),
           ZNF501 = round(`ZNF501|115560`)) -> ACC_BLCA_BRCA_OV.rnaseq.medians

# 以癌症种类为列，以"MET|4233"表达区间为行，以"ZNF500|26048"的区间均值为值，绘图
heatmapTCGA(ACC_BLCA_BRCA_OV.rnaseq.medians, "cohort", "MET", "ZNF500", title = "Heatmap of ZNF500 expression")
