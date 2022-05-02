a = read.csv("r-pkg.xls", sep="\t")

list.of.packages <- a[,1]

list.of.bio.packages  <- a[,1]

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos="http://cran.rstudio.com/", dependencies=TRUE)

new.bio.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.bio.packages)){
    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install(new.bio.packages, suppressUpdates = T)
}	

WGCNA must be installed after AnnotationDbi, impute, GO.db, preprocessCore
install.packages(c("WGCNA"))

Test Load Packages    This cause an error when too many packages are loaded.
suc = unlist ( lapply(list.of.packages, require, character.only = TRUE) )
if(sum(suc) < length(list.of.packages) )
	cat ("\n\nWarnning!!!!!! These R packages cannot be loaded:", list.of.packages[!suc] )

suc = unlist ( lapply(list.of.bio.packages, require, character.only = TRUE) )
if(sum(suc) < length(list.of.bio.packages) )
	cat ("\n\nWarnning!!!!!! These Bioconductor packages cannot be loaded:", list.of.bio.packages[!suc] )