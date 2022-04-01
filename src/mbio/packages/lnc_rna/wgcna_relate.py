# coding=utf-8
import subprocess
import argparse

"""
1. Correlate the module eigengenes with the trait
2. Correlate each gene with the trait
"""
# get arguments
parser = argparse.ArgumentParser(description="wgcna step3: relate module/gene to traits")
parser.add_argument('-datExpr', type=str, metavar="exp_matrix_file", required=True,
                    help="expression matrix file, for gene/traits relationship calculation")
parser.add_argument('-MEs', type=str, metavar="module eigengenes", required=True,
                    help="module eigengenes, for module/traits relationship calculation")
parser.add_argument('-traits', type=str, required=True, metavar="phenotype_data",
                    help="sample name in row, traits data in column")
parser.add_argument('-corType', type=str, metavar="correlation_type", default="pearson",
                    help="correlation type, 'pearson' or 'spearman', 'kendall'")
parser.add_argument('-nThreads', type=int, default=16, )
parser.add_argument('-block_Rdata', type=str, default=None, )
args = parser.parse_args()

# read expr
r_cmds = \
"""
# module and traits relation
library('WGCNA')
enableWGCNAThreads()
datME = read.table('{eigengenes}', header=T, sep="\\t", row.names=1)
traits = read.table("{traits}", header=T, sep="\\t", row.names=1)
if (dim(traits)[2] == 1 & class(traits[1,1])=="factor"){bracket1} 
    tmp = model.matrix(~0+ traits[,1])
    colnames(tmp) = levels(traits[,1])
    traits = tmp
{bracket2}
correlation = signif(cor(t(datME), traits, use="p", method="{cor_type}", nThreads={threads}), 3)
pvalues = signif(corPvalueStudent(correlation, nSamples = dim(traits)[1]), 3)
write.table(correlation, 'module_trait.correlation.xls', col.names=NA, quote=F, sep='\\t', row.names=T)
write.table(pvalues, 'module_trait.correlation_pvalues.xls', col.names=NA, quote=F, sep='\\t', row.names=T)
# gene and traits relation
exp = read.table('{exp_matrix}', header=T, row.names=1)
correlation2 = signif(cor(t(exp), traits, use="p", method="{cor_type}"), 3)
write.table(correlation2, 'gene_trait.correlation.xls', col.names=NA, quote=F, sep='\\t', row.names=T)
traitColors = numbers2colors(correlation2, signed = FALSE, colors=greenWhiteRed(100))

load("{rdata}")
block_number = length(bwnet$dendrograms)
for (i in c(1:block_number)){bracket1}
block_traitColors = traitColors[bwnet$blocks==i,]
all_color = cbind(bwnet$color[bwnet$blockGenes[[i]]], block_traitColors)
labels = c(c("Module colors"), colnames(correlation2))
pdf(paste('block_', i, '_gene_trait.correlation.pdf', sep=''), width=12, height=9)
plotDendroAndColors(bwnet$dendrograms[[i]], all_color, labels, main = paste("Gene dendrogram and module colors in block", 1),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
png(paste('block_', i, '_gene_trait.correlation.png', sep=''), width=1200, height=900, type = "cairo")
plotDendroAndColors(bwnet$dendrograms[[i]], all_color, labels, main = paste("Gene dendrogram and module colors in block", 1),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
{bracket2}
""".format(
    eigengenes=args.MEs,
    traits=args.traits,
    cor_type=args.corType,
    threads=args.nThreads,
    exp_matrix=args.datExpr,
    bracket1="{",
    bracket2="}",
    rdata=args.block_Rdata
)

with open('wgcna_relate_analysis.r', 'w') as f:
    f.write(r_cmds)
subprocess.check_call("Rscript wgcna_relate_analysis.r", shell=True)
