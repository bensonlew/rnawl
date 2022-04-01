# Load the bioconductor installer.
source("https://bioconductor.org/biocLite.R")
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

# Install the main RTCGA package
biocLite("RTCGA")

# Install the clinical and mRNA gene expression data packages
biocLite("RTCGA.clinical") ## (14.0 MB)
biocLite("RTCGA.rnaseq") ## (612.6 MB)
biocLite("RTCGA.mRNA") ## (85.0 MB)
biocLite("RTCGA.mutations") ## (103.8 MB)

library(RTCGA)
library(RTCGA.clinical)
library(RTCGA.rnaseq)
library(RTCGA.mRNA)
library(RTCGA.mutations)

all_TCGA_cancers <- infoTCGA()
