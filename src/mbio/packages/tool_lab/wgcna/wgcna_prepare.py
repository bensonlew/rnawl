import json
import argparse
import pandas as pd
import numpy as np
import subprocess
import fastcluster as hclust
import scipy.cluster.hierarchy as sch


# get arguments
parser = argparse.ArgumentParser(description="wgcna step1: pre-processing data")
parser.add_argument('-exp', type=str, metavar="exp_matrix_file", required=True,
                    help="expression matrix file")
parser.add_argument('-me', type=float, metavar="mean_exp_cutoff", default=None,
                    help="mean expression cutoff")
parser.add_argument('-cv', type=float, metavar="coefficient_variance_cutoff", default=None,
                    help="mean expression cutoff")

args = parser.parse_args()
exp_file = args.exp
low_exp_cutoff = args.me
cv_cutoff = args.cv

# clean data
try:
    exp_pd = pd.read_table(exp_file, header=0, index_col=0)
except:
    exit('pandas failed')
exp_pd.index = exp_pd.index.astype('str', copy=False)   # added by zhangyitong on 20210816
exp_pd = np.log2(exp_pd+1)
initial_genes = exp_pd.index
if low_exp_cutoff:
    exp_pd = exp_pd[exp_pd.mean(axis=1) > low_exp_cutoff]
if cv_cutoff:
    exp_pd = exp_pd[exp_pd.std(axis=1)/exp_pd.mean(axis=1) > cv_cutoff]
if exp_pd.shape[0] == 0:
    exit("oh, NO! All genes are filtered")
final_genes = exp_pd.index
ignored_genes = set(initial_genes) - set(final_genes)
with open('ignored_gene.list', 'w') as f:
    f.write('# {}/{}\n'.format(len(final_genes), len(initial_genes)))
    for each in ignored_genes:
        f.write(each + '\n')
exp_pd.to_csv('exp_matrix_after_filtering.txt', sep='\t', index=True, header=True)

# cluster sample using python
# z = hclust.linkage(exp_pd.transpose(), method="average", metric="euclidean")
# min_h = min(z[:, 2])
# dendrogram = sch.dendrogram(z, labels=exp_pd.columns, distance_sort='ascending')
# with open('dendrogram.json', 'w') as f:
#     json.dump(dendrogram, f)
# with open('dendrogram.json', 'a+') as f:
#     f.write("\n{}\n".format(min_h))


# pick power and cluster sample
r_cmd1 = [
    "library('WGCNA')",
    "enableWGCNAThreads()",
    "exp_matrix = read.table('exp_matrix_after_filtering.txt', header=T, row.names=1, sep='\\t')",
    "datExpr = as.data.frame(exp_matrix)",
    "datExpr = t(datExpr)",
]

# cluster sample using R
r_cmd1.append(
"""
sampleTree = hclust(dist(datExpr), method = "average");
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
linkage = cbind(sampleTree$merge, sampleTree$height)
write.table(linkage, 'sample.cluster.dendrogram.txt', col.names=F, quote=F, sep='\\t', row.names=F)
write.table(sampleTree$order, 'sample.cluster.dendrogram.order.txt', col.names=F, quote=F, sep='\\t', row.names=F)
write.table(sampleTree$labels[sampleTree$order], 'sample.cluster.dendrogram.labels.txt', col.names=F, quote=F, sep='\\t', row.names=F)
""")

r_cmd2 = ["""
# calculate powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
target_power = sft$powerEstimate
cutoff = sft$fitIndices[target_power, 2]
# save result
write.table(sft$fitIndices, 'scale_free_analysis.xls', sep='\t', quote=F, col.names=T, row.names=F)
write.table(cutoff, paste('powerEstimate_', target_power, sep=''), row.names=F, col.names=F)

# Plot the results:
pdf(file='pick_power.pdf', width = 9, height = 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=cutoff, col="green")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
"""]

final_cmd = "\n".join(r_cmd1 + r_cmd2)
with open('wgcna.r', 'w') as f:
    f.write(final_cmd)
subprocess.check_call("Rscript wgcna.r", shell=True)

