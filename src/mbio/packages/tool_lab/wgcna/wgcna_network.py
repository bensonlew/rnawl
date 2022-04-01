import subprocess
import argparse
import glob
import json
import os
import pandas as pd
import numpy as np


# get arguments
parser = argparse.ArgumentParser(description="wgcna step4: get network data")
parser.add_argument('-module', type=str, metavar="target_module", required=True, help="such as: blue,red")
parser.add_argument('-threshold', type=float, metavar="weight_threshold", default=0.02)
parser.add_argument('-top', type=int, metavar="nTop_hubGenes", default=30)
parser.add_argument('-step2output', type=str, metavar="output_dir", required=True,
                    help="NeedFiles: gene_trait.correlation.xls")
parser.add_argument('-step3output', type=str, metavar="output_dir", required=True,
                    help="NeedFiles: TOM-blockwise-block.*.RData, seq_id2gene_name.txt, "
                         "blockwiseModules_result.RData, membership.xls")
args = parser.parse_args()

tom_files = glob.glob(args.step2output+"/TOM-blockwise-block.*.RData")
bwm_result = glob.glob(args.step2output+"/blockwiseModules_result.RData")[0]
gene_annot = glob.glob(args.step2output+"/seq_id2gene_name.txt")[0]
if not tom_files:
    raise Exception("No TOM files found")
target_modules = args.module.strip().split(",")

# format R scripts
r_cmd = ["""
library("WGCNA")
load('{bwnet}')
moduleColors = membership[, 2];
probes = membership[, 1];
annot = read.table('{gene_annot}', header=T, sep='\\t');
""".format(gene_annot=gene_annot, bwnet=bwm_result, )
]

for tom_file in tom_files:
    block_id = os.path.basename(tom_file).split(".")[-2]
    r_cmd.append(
"""
# ---export a block network---
load('{tom_file}');
blockTOM = as.matrix(TOM);
blockModuleColors = moduleColors[bwnet$blockGenes[[{block_id}]]];
inModule = is.finite(match(blockModuleColors, c({modules})));
if (length(which(inModule == TRUE)) >= 4){bracket1}
    blockProbes = probes[bwnet$blockGenes[[{block_id}]]];
    modTOM = blockTOM[inModule, inModule];
    modProbes = blockProbes[inModule];
    dimnames(modTOM) = list(modProbes, modProbes);
    modGenes = annot$gene_name[match(modProbes, annot[,1])];
    modColors = blockModuleColors[inModule];
    nodeAttrs = modColors;
    # get top n hub genes
    IMConn = softConnectivity(datExpr[, modProbes]);
    top = (rank(-IMConn, ties.method='first') <= {top})
    modTOM = modTOM[top, top]
    modProbes = modProbes[top]
    modGenes = modGenes[top]
    modColors = modColors[top]
    nodeAttrs = cbind(modColors, IMConn[top])
    colnames(nodeAttrs) = c('module', 'intraModuleConnectivity')
    # to cytoscape
    cyt = exportNetworkToCytoscape(modTOM,
      edgeFile = "b{block_id}_edges.txt",
      nodeFile = "b{block_id}_nodes.txt",
      weighted = TRUE,
      threshold = {threshold},
      nodeNames = modProbes,
      altNodeNames = modGenes,
      nodeAttr = nodeAttrs);
{bracket2}
""".format(
    tom_file=tom_file,
    block_id=block_id,
    modules=",".join("'"+x+"'" for x in target_modules),
    threshold=args.threshold,
    top=args.top,
    bracket1="{",
    bracket2="}"
))
else:
    r_file = 'wgcna_network.r'
    with open(r_file, 'w') as f:
        f.write('\n'.join(r_cmd))
    subprocess.check_call("Rscript {}".format(r_file), shell=True)

# merge network
all_edges = glob.glob("b*_edges.txt")
all_nodes = glob.glob("b*_nodes.txt")
# edge
links = list()
for each in all_edges:
    edge_pd = pd.read_table(each, header=0)
    # "fromNode	toNode	weight	direction	fromAltName	toAltName"
    if edge_pd.shape[0] <= 2:
        continue
    edge_pd.columns = ['source', 'target', 'weight', 'direction', 'fromAltName', 'toAltName']
    links.extend(edge_pd.to_dict("records"))

if links:
    # exit("Network is Empty! Stop!")
    # node
    gene_trait_corr = glob.glob(args.step3output+"/gene_trait.correlation.xls")[0]
    gene_trait_corr_pd = pd.read_table(gene_trait_corr, header=0, index_col=0)
    gene_trait_corr_pd.columns = [x + "_corr" for x in gene_trait_corr_pd.columns]
    gene_trait_corr_pd.index.name = "seq_id"
    gene_module_corr = glob.glob(args.step2output+"/membership.xls")[0]
    gene_module_corr_pd = pd.read_table(gene_module_corr, header=0, index_col=0).loc[:, ["kme"]]
    gene_module_corr_pd.index.name = "seq_id"
    gene_annot_pd = pd.read_table(gene_annot, header=0, index_col=0)
    if 'gene_id' in gene_annot_pd.columns:
        gene_annot_pd = gene_annot_pd.loc[:, ["gene_id"]]
    else:
        gene_annot_pd.index.name = 'seq_id'
        gene_annot_pd['gene_id'] = gene_annot_pd.index
        if "gene_type" in set(gene_annot_pd.columns):
            gene_annot_pd = gene_annot_pd.loc[:, ["gene_id", "gene_type"]]
        else:
            gene_annot_pd = gene_annot_pd.loc[:, ["gene_id"]]
    nodes = list()
    for each in all_nodes:
        node_pd = pd.read_table(each, header=0, index_col=0)
        node_pd.columns = ['name', 'module', 'connectivity']
        node_pd.index.name = "seq_id"
        node_pd = node_pd.join(gene_trait_corr_pd.round(4), how="left")
        node_pd = node_pd.join(gene_module_corr_pd.round(4), how="left")
        if 'gene_id' not in list(node_pd.columns):
            node_pd = node_pd.join(gene_annot_pd, how="left")
        node_pd.index.name = "id"
        node_pd.reset_index(inplace=True)
        # nodeName	altName	nodeAttr.nodesPresent...
        nodes.extend(node_pd.to_dict("records"))

    if len(all_nodes) > 1:
        # multi block may have multi same modude
        nodes = sorted(nodes, key=lambda x:x['connectivity'], reverse=True)[:args.top]
        node_list = [x['id'] for x in nodes]
        links = [x for x in links if x['source'] in node_list and x['target'] in node_list]
    # write out network for cytoscape
    link_order = ['source', 'target', 'weight', 'direction', 'fromAltName', 'toAltName']
    final_link_pd = pd.DataFrame(links).loc[:, link_order]
    final_link_pd = final_link_pd.fillna("-")
    # final_link_pd.to_csv("_".join(target_modules)+".network.edges.txt", header=True, index=False, sep='\t')
    final_link_pd.to_csv("network.edges.txt", header=True, index=False, sep='\t')
    final_node_pd = pd.DataFrame(nodes)
    node_order = ["id", "name", "connectivity", "gene_id", "kme", "module"] + list(set(final_node_pd.columns) - {"id","gene_id", "name","kme","module","connectivity"})
    final_node_pd.loc[final_node_pd['name'].isnull(), 'name'] = final_node_pd[final_node_pd['name'].isnull()]['id']
    # final_node_pd.loc[:, node_order].to_csv("_".join(target_modules)+".network.nodes.txt", header=True, index=False, sep='\t')
    final_node_pd.loc[:, node_order].to_csv("network.nodes.txt", header=True, index=False, sep='\t')
    # write out network for mj sanger
    # with open("_".join(target_modules)+".network.json", "w") as f:
    with open("network.json", "w") as f:
        tmp_nodes = final_node_pd.loc[:, ["id", "name", "module"]]
        tmp_nodes.columns = ["id", "name", "group"]
        index_node_pd = tmp_nodes.loc[:, ["id"]].reset_index().set_index("id")
        tmp_links = final_link_pd.loc[:, ["source", "target", "weight"]].round(4)
        tmp_links.columns = ["source", "target", "distance"]
        tmp_links = tmp_links.join(index_node_pd, on="source")
        tmp_links = tmp_links.join(index_node_pd, on="target", lsuffix="_source", rsuffix="_target")
        tmp_links = tmp_links.loc[:, ["index_source", "index_target", "distance"]]
        tmp_links.columns = ["source", "target", "distance"]
        json.dump(dict(links=json.loads(tmp_links.to_json(orient="records")), nodes=tmp_nodes.to_dict("records")), f)
