# load package

library(getopt)

# set options
command <- matrix(c(
    "json", "i", 1, "character", "setting file in JSON format",
    "help", "h", 0, "logical", "show this help message and exit"
), byrow = TRUE, ncol = 5)
opts <- getopt(command)

if (! is.null(opts$help)) {
    cat(getopt(command, usage = TRUE))
    q(status = 2)
}
if (is.null(opts$json)) {
    cat(getopt(command, usage = TRUE))
    q(status = - 1)
}

library(rjson)
library(edgeR)
library(limma)
#library(DESeq2)


# set variables

print("start loading required data")
json.lst <- fromJSON(file = opts$json)
count.path <- json.lst$count
group.vec <- json.lst$groups
BATCH = json.lst$Batch
PAIRED = json.lst$Paired
batch.vec = json.lst$batchs
paired.vec = json.lst$paireds
contrast.vec <- json.lst$contrasts
output.path <- json.lst$output
#fit.method <- json.lst$fit_method
norm = json.lst$norm_method #Quantile_Normalization,calcNormFactors,calcNormFactors+voom,voomWithQualityWeights,voom+vst
trend = json.lst$trend # TRUE or FALSE
robust = json.lst$robust #TRUE or FALSE
count.df <- read.delim(count.path, row.names = 1)
count.df = floor(count.df)

print("start creating group.factor")
group = as.factor(group.vec)

if (BATCH == TRUE) {
print("start creating batch.factor")
batch = as.factor(batch.vec)
}

if (PAIRED == TRUE) {

print("start creating paired.factor")
paired = as.factor(paired.vec)
}

print("start creating design matrix")
if (BATCH == FALSE & PAIRED == FALSE) {
design.mat = model.matrix(~0+group)
design.mat
#colnames(design.mat)=levels(group)
#rownames(design.mat) <- colnames(count.df)
} else if (BATCH == TRUE & PAIRED == FALSE) {
design.mat = model.matrix(~0+group+batch)
design.mat
#colnames(design.mat)=levels(group)
#rownames(design.mat) <- colnames(count.df)
} else if (BATCH == FALSE & PAIRED == TRUE) {
design.mat = model.matrix(~0+group+paired)
design.mat
}

print("start Computational differential expression")
if (norm == 'quantile') {
  counts.log=log2(count.df+1)
  counts.log.norm=normalizeBetweenArrays(counts.log,method='quantile')
  fit = lmFit(counts.log.norm,design.mat)

} else if (norm == 'calcNormFactors') {
    dgel <- DGEList(counts=count.df,group=group)
    dgel <- calcNormFactors(dgel)
    logCPM <- cpm(dgel, log=TRUE)
    fit = lmFit(logCPM,design.mat)
} else if (norm == 'cal_voom') {
    dgel = DGEList(counts=count.df,group=group)
    dgel = calcNormFactors(dgel)
    dgel_voom = voom(dgel,design.mat)
    fit = lmFit(dgel_voom,design.mat)
} else if (norm == 'cal_voomqualityweight') {
    dgel = DGEList(counts=count.df,group=group)
    dgel = calcNormFactors(dgel)
    dgel.voom = voomWithQualityWeights(dgel,design.mat,normalization="none",plot=FALSE)
    fit = lmFit(dgel.voom,design.mat)
} else if (norm == 'voom_vst') {
    cds <- newCountDataSet(count.df, group)
    cds <- estimateSizeFactors(cds)
    cds <- estimateDispersions(cds,method='blind',fitType='local')
    vst <- getVarianceStabilizedData(cds)
    fit <- lmFit(vst,design.mat)
}

for (i in seq(length(contrast.vec))) {
  case = contrast.vec[[i]]$case
  ctrl = contrast.vec[[i]]$ctrl
  group_case = paste('group',case,sep = "")
  group_ctrl = paste('group',ctrl,sep = "")
  cont.matrix = makeContrasts(paste(group_case, group_ctrl, sep = "-"),levels = design.mat)
  fit2 = contrasts.fit(fit,cont.matrix)
  fit3 = eBayes(fit2,trend = trend,robust = robust)
  results = topTable(fit3,n=nrow(count.df))
  o.path <- file.path(output.path, paste(paste(ctrl, case, sep = "_vs_"), "limma", "txt", sep = "."))
  write.table(t(c("seq_id", colnames(results))), file = o.path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(results, file = o.path, append = TRUE, quote = FALSE, sep = "\t", col.names = FALSE)
  print(paste("succeed in exporting", o.path))

  }

