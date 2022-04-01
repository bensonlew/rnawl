# BiocManager::install("KEGG.db")
library(cgdsr)
library(KEGG.db)

# 获取细胞循环和细胞因子KEGG通路内基因的entrez编号向量
cellcycle.gene.entrez <- KEGGPATHID2EXTID[['hsa04110']]
cytokine.gene.entrez <- KEGGPATHID2EXTID[['hsa04060']]

# 创建CGDS连接对象
mycgds <- CGDS("http://www.cbioportal.org/public-portal/")

cancer.studies <- getCancerStudies(mycgds)