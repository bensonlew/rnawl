# load package
library(Rsubread)

args = commandArgs(T)
print(args)
bam <- args[1]
gtf <- args[2]
isPairedEnd <- args[3]
countMultiMappingReads <- args[4]
allowMultiOverlap <-args[5]
fraction <- args[6]
output <- args[7]
fc <- featureCounts(bam, annot.ext=gtf, isGTFAnnotationFile=TRUE, GTF.featureType='exon', GTF.attrType="gene_id", useMetaFeatures=TRUE, isPairedEnd=isPairedEnd, countMultiMappingReads=countMultiMappingReads, allowMultiOverlap=allowMultiOverlap, fraction=fraction, strandSpecific=1)
counts <- fc$counts
write.table(counts, file=output, sep='\t', col.names=FALSE, quote=FALSE)







