# coding=utf-8
import subprocess
import argparse
import pandas as pd
import glob
import os

# get arguments
parser = argparse.ArgumentParser(description="wgcna step1: pre-processing data")
parser.add_argument('-datExpr', type=str, metavar="exp_matrix_file", required=True,
                    help="expression matrix file")
parser.add_argument('-mergeCutHeight', type=float, metavar="merge_cut_height", default=0.25,
                    help="dynamic_cut_tree height cutoff")
parser.add_argument('-corType', type=str, metavar="correlation_type", default="pearson",
                    help="correlation type, 'pearson' or 'bicor'")
parser.add_argument('-maxBlockSize', type=int, default=20000, metavar="max_block_size",
                    help="maximum block size for block wise module")
parser.add_argument('-power', type=int, default=6,
                    help="power value picked for Adjacency function options")
parser.add_argument('-networkType', type=str, default='signed', metavar="network_type",
                    help="network type, signed or unsigned or signed hybrid")
# parser.add_argument('-TOMType', type=str, default='signed', metavar="TOM_type",
#                     help="TOM type, signed or unsigned or none")
parser.add_argument('-minModuleSize', type=int, default=0, metavar="min_module_size",
                    help="min module size for block wise module. default: min(20, ncol(datExpr)/2 )")
parser.add_argument('-minKMEtoStay', type=float, default=0.3, )
parser.add_argument('-nThreads', type=int, default=16, )
args = parser.parse_args()

# read expr
exp_matrix = pd.read_table(args.datExpr, header=0, index_col=0)
# special args
if args.minModuleSize == 0:
    gene_number = exp_matrix.shape[0]
    min_module = min([20, int(gene_number/2)])
else:
    min_module = args.minModuleSize
if "unsigned" in args.networkType:
    tom_type = "unsigned"
else:
    tom_type = "signed"

# step1: read data
r_cmd1 = [
    "library('WGCNA')",
    "library('ape')",
    "enableWGCNAThreads()",
    "exp_matrix = read.table('{exp}', check.names = FALSE, header=T, row.names=1, sep='\\t')".format(exp=args.datExpr),
    "datExpr = as.data.frame(exp_matrix)",
    "datExpr = t(datExpr)",
]

# step2: blockwiseModules
r_cmd2 = ['bwnet = blockwiseModules('
          'datExpr, '
          'maxBlockSize = {block_size}, '
          'power = {power}, '
          'minModuleSize = {min_module},'
          'mergeCutHeight = {cut_height}, '
          'nThreads={threads}, '
          'numericLabels = FALSE, '
          'saveTOMs = TRUE, '
          'corType = "{corr_type}",'
          'saveTOMFileBase = "TOM-blockwise", '
          'networkType = "{net_type}", '
          'TOMType = "{tom_type}", '
          'minKMEtoStay = {kme_stay}, '
          'verbose = 3)'.format(block_size=args.maxBlockSize,
                                power=args.power,
                                threads=args.nThreads,
                                min_module=min_module,
                                cut_height=args.mergeCutHeight,
                                corr_type=args.corType,
                                net_type=args.networkType,
                                tom_type=tom_type,
                                kme_stay=args.minKMEtoStay,
                                )
]

# step3: Plot the dendrogram and the module colors underneath
r_cmd3 = ["""
for (i in c(1:length(bwnet$dendrograms))){
pdf(file=paste('block_', i, '_dendrogram.pdf', sep=''), width = 15, height = 9)
plotDendroAndColors(bwnet$dendrograms[[i]],bwnet$color[bwnet$blockGenes[[i]]],"Module colors",
main = paste("Gene dendrogram and module colors in block", i),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
png(file=paste('block_', i, '_dendrogram.png', sep=''), width=1080, height=390, type = "cairo")
plotDendroAndColors(bwnet$dendrograms[[i]],bwnet$color[bwnet$blockGenes[[i]]],"Module colors",
main = paste("Gene dendrogram and module colors in block", i),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
}
"""]

# export module and block info
r_cmd4 = ["""
# save color list
membership = cbind(colnames(datExpr), bwnet$colors, bwnet$blocks)
write.table(membership, 'color.list', col.names=F, quote=F, sep='\\t', row.names=F)
# save dendrogram for each block
block_number = length(bwnet$dendrograms)
for (i in c(1:block_number)){
    write.table(t(bwnet$MEs), 'eigengenes.txt', sep='\\t', quote=F, col.names=row.names(datExpr),)
    block_dendrogram = bwnet$dendrograms[[i]]
    linkage = cbind(block_dendrogram$merge, block_dendrogram$height)
    write.table(linkage, paste('block_', i, '.dendrogram.txt', sep=''), col.names=F, quote=F, sep='\\t', row.names=F)
    block_dendrogram_order = block_dendrogram$order
    write.table(block_dendrogram_order, paste('block_', i, '.dendrogram.order.txt', sep=''), col.names=F, quote=F, sep='\\t', row.names=F)
    labels = bwnet$color[bwnet$blockGenes[[i]][block_dendrogram_order]]
    ordered_genes = colnames(datExpr)[bwnet$blockGenes[[i]][block_dendrogram_order]]
    write.table(labels, paste('block_', i, '.dendrogram.labels.txt', sep=''), col.names=F, quote=F, sep='\\t', row.names=F)
    write.table(ordered_genes, paste('block_', i, '.dendrogram.ordered_seqs.txt', sep=''), col.names=F, quote=F, sep='\\t', row.names=F)
}
"""]

# get membership
r_cmd5 = ["""
kme = signedKME(datExpr, bwnet$MEs)
write.table(kme, 'kme.table.xls', col.names=NA, quote=F, sep='\\t', row.names=T)
"""]

# 模块相关性聚类分析
r_cmd6 = ['MEs = bwnet$MEs', 'MECorr = cor(MEs, method="p")']
if args.networkType == 'signed':
    r_cmd6.append('adjacency = (1+MECorr)/2')
    r_cmd6.append('MEDiss = (1-MECorr)/2')
else:
    r_cmd6.append('adjacency = abs(MECorr)')
    r_cmd6.append('MEDiss = 1-abs(MECorr)')
r_cmd6.append("""
METree = hclust(as.dist(MEDiss), method = "average")
tree = as.phylo(METree)
write.tree(tree, 'module_corr.tree.txt')
write.table(adjacency, 'module_corr.matrix.xls', col.names=NA, quote=F, sep='\\t', row.names=T)
""")

# 针对模块合并前的eigengens计算聚类树
r_cmd7 = ["""
MEs0 <- moduleEigengenes(datExpr, bwnet$unmergedColors)$eigengenes
MECorr = cor(MEs0, method="p")
"""]
# if args.networkType == 'signed':
#     r_cmd7.append('adjacency = (1+MECorr)/2')
#     r_cmd7.append('MEDiss = (1-MECorr)/2')
# else:
#     r_cmd7.append('adjacency = abs(MECorr)')
#     r_cmd7.append('MEDiss = 1-abs(MECorr)')
r_cmd7.append('adjacency = abs(MECorr)')
r_cmd7.append('MEDiss = 1-abs(MECorr)')
r_cmd7.append("""
METree = hclust(as.dist(MEDiss), method = "average")
error_info = try({bracket1}
pdf(file="cluster_unmerged_module_eigengenes.pdf", width = 17, height = 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "", cex=0.6)
abline(h={height}, col = "red")
dev.off()
{bracket2}, silent=T)
print(error_info)
linkage = cbind(METree$merge, METree$height)
write.table(linkage, 'unmerged.module_corr.dendrogram.txt', col.names=F, quote=F, sep='\\t', row.names=F)
write.table(METree$order, 'unmerged.module_corr.dendrogram.order.txt', col.names=F, quote=F, sep='\\t', row.names=F)
write.table(METree$labels[METree$order], 'unmerged.module_corr.dendrogram.labels.txt', col.names=F, quote=F, sep='\\t', row.names=F)
save(bwnet, membership, datExpr, file = "blockwiseModules_result.RData")
""".format(height=args.mergeCutHeight, bracket1="{", bracket2="}"))






final_cmd = "\n".join(r_cmd1 + r_cmd2 + r_cmd3 + r_cmd4 + r_cmd5 + r_cmd6 + r_cmd7)
with open('wgcna.r', 'w') as f:
    f.write(final_cmd)
subprocess.check_call("Rscript wgcna.r", shell=True)

tom_files = glob.glob("TOM-blockwise-block.*.RData")
bwm_result = glob.glob("blockwiseModules_result.RData")[0]
r_cmd = ["""
library("WGCNA")
load('{bwnet}')
moduleColors = membership[, 2];
probes = membership[, 1];
""".format(bwnet=bwm_result, )
]
for tom_file in tom_files:
    block_id = os.path.basename(tom_file).split(".")[-2]
    r_cmd.append(
"""
# ---export a block network---
load('{tom_file}');
blockTOM = as.matrix(TOM);
blockModuleColors = moduleColors[bwnet$blockGenes[[{block_id}]]];
dissTOM=1-blockTOM;
selectTOM = dissTOM;
selectTree = hclust(as.dist(selectTOM),method = "average");
selectColors = moduleColors;
plotDiss=selectTOM^7;
diag(plotDiss)=NA;
pdf_name = '{block_id}_Network-heatmap.pdf';
pdf(pdf_name, width = 15, height = 9);
TOMplot(plotDiss,selectTree,selectColors,main="Network heapmap");
dev.off();
png_name = '{block_id}_Network-heatmap.png';
png(png_name, width=1080, height=600, type = "cairo")
TOMplot(plotDiss,selectTree,selectColors,main="Network heapmap");
dev.off();
""".format(
    tom_file=tom_file,
    block_id=block_id,
))

with open('wgcna_1.r', 'w') as f:
    f.write('\n'.join(r_cmd))
subprocess.check_call("Rscript wgcna_1.r", shell=True)

color_list = pd.read_table('color.list', header=None, dtype={0: 'str'})
gene2module = dict(zip(color_list[0], color_list[1]))
gene2block = dict(zip(color_list[0], color_list[2]))
kme_pd = pd.read_table('kme.table.xls', header=0, index_col=0)
kme_pd.index = kme_pd.index.astype('str', copy=False)   # added by zhangyitong on 20210816
gene2kme = dict()
for each in kme_pd.index:
    if gene2module.has_key(each):
        gene2kme[each] = kme_pd.loc[each, 'kME'+gene2module[each]]
    elif "." in each:
        each_new = each.replace(".", "-")
        if gene2module.has_key(each_new):
            gene2kme[each_new] = kme_pd.loc[each, 'kME'+gene2module[each_new]]
    else:
        raise Exception("cannot find {} in gene list".format(each))
with open('membership.xls', 'w') as f:
    f.write("{}\t{}\t{}\t{}\n".format('seq_id', 'module', 'kme', 'block_id'))
    for each in color_list[0]:
        f.write('{}\t{}\t{}\t{}\n'.format(each, gene2module[each], gene2kme[each], gene2block[each]))
# module stat
groups = color_list.groupby(1).groups
with open('module_size.stat.xls', 'w') as f:
    f.write("module\tsize\tseqs\n")
    for each in groups:
        genes = color_list.iloc[groups[each], :][0]
        f.write("{}\t{}\t{}\n".format(each, len(genes), ";".join(genes)))


# ------------------plot-----------------
import pandas as pd
import glob
import matplotlib.pyplot as plt


def get_tree_coord(linkage_pd, ordered_leaf, ):
    """
    format block tree data
    :param ordered_leaf: [1,2,3,..],代表画聚类树的结果
    :param linkage_pd: [n1,n2,height] 来自R的聚类结果，不适用python的聚类结果
    :return:
    """
    sn = len(ordered_leaf)
    number = sn
    leaf_height_dict = {x: 0 for x in range(sn)}
    origin_leaf_x_loc = range(5, (sn - 1) * 10 + 6, 10)
    leaf_x_dict = dict(zip(ordered_leaf, origin_leaf_x_loc))
    x_coord = list()
    y_coord = list()
    steps = list()
    min_height = linkage_pd.iloc[0, 2]
    for each in linkage_pd.iterrows():
        sn += 1
        n1, n2 = int(each[1][0]), int(each[1][1])
        if n1 > 0:
            n1 = number + n1
        else:
            n1 = abs(n1)
        if n2 > 0:
            n2 = number + n2
        else:
            n2 = abs(n2)
        dis = float(each[1][2])
        leaf_height_dict[sn] = dis
        leaf_x_dict[sn] = (leaf_x_dict[n1] + leaf_x_dict[n2]) / 2
        if n1 <= number:
            y1 = dis - min_height*0.01
        else:
            y1 = leaf_height_dict[n1]
        y2, y3 = dis, dis
        if n2 <= number:
            y4 = dis - min_height*0.01
        else:
            y4 = leaf_height_dict[n2]
        #     y1 = leaf_height_dict[n1]
        #     y4 = leaf_height_dict[n2]
        steps.append([n1, n2])
        y_coord.append([y1, y2, y3, y4])
        x_coord.append([leaf_x_dict[n1], leaf_x_dict[n1], leaf_x_dict[n2], leaf_x_dict[n2]])
    return x_coord, y_coord, ordered_leaf


def plot_gene_tree(result_dir):
    blocks = glob.glob(result_dir + '/block_*.dendrogram.txt')
    fig, axes = plt.subplots(len(blocks), 1)
    for each in blocks:
        linkage_pd = pd.read_table(each, header=None)
        ordered_leaf = pd.read_table(each[:-3] + 'order.txt', header=None)[0]
        x_coord, y_coord,  _ = get_tree_coord(linkage_pd, ordered_leaf)
        for x, y in zip(x_coord, y_coord):
            if len(blocks) >= 2:
                axes[0].plot(x, y)
            else:
                axes.plot(x, y)
    plt.savefig('gene_tree.png')


# plot_gene_tree(".")
