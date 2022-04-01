#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import pandas as pd
import json
import time
import os

parser = argparse.ArgumentParser()
parser.add_argument('-input_file', type=str, required=True,help='path of input file')
parser.add_argument('-cosmic_file_path', type=str, required=True,help="cosmic_file_path")
parser.add_argument('-output', type=str, required=True,help="output file path")
parser.add_argument('-pub_json', type=str, required=True,help="pubmed json path")
args = parser.parse_args()
output = args.output
input_file = args.input_file
cosmic_file_path = args.cosmic_file_path
pub_json = args.pub_json

fusions = []
cosmic_df = pd.read_table(cosmic_file_path)
with open(input_file,"r") as r:
    r.readline()
    for line in r.readlines():
        left_gene = line.strip().split("\t")[0]
        right_gene = line.strip().split("\t")[1]
        fusions.append((left_gene,right_gene))

target_dfs = []
for fusion in fusions:
    left_gene,right_gene = fusion
    target_df = cosmic_df[(cosmic_df["5'_GENE_NAME"].str.contains(left_gene)) & (cosmic_df["3'_GENE_NAME"].str.contains(right_gene))]
    target_df["Source"] = "Cosmic"
    target_df["Gene"] = left_gene + "_" + right_gene
    target_df["Genome"] = "hg19"
    target_df_reverse = cosmic_df[(cosmic_df["3'_GENE_NAME"].str.contains(left_gene)) & (cosmic_df["5'_GENE_NAME"].str.contains(right_gene))]
    target_df_reverse["Source"] = "Cosmic"
    target_df_reverse["Gene"] = left_gene + "_" + right_gene
    target_df_reverse["Genome"] = "hg19"
    target_dfs.append(target_df)
    target_dfs.append(target_df_reverse)
if target_dfs:
    final_df = pd.concat(target_dfs)
    final_df = final_df[["Source","Gene","Genome","PRIMARY_HISTOLOGY","PUBMED_PMID"]]
    pub2name = json.load(open(pub_json))
    final_df[" references_title"] = final_df["PUBMED_PMID"].apply(lambda x : pub2name[str(x)] if str(x) in pub2name else "-")
    final_df.columns = ["Source","Gene","Genome","Disease","PUBMED_PMID","Reference"]
    final_df.to_csv(os.path.join(output,"Fusion_anno.xls"),sep = "\t",index =False)
else:
    with open(os.path.join(output,"Fusion_anno.xls"),"w") as w:
        w.write("无注释结果")



