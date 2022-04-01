# -*- coding: utf-8 -*-
import pandas as pd
from biocluster.config import Config
import os
import argparse

def select_genes(statfile, outfile, select_sam=None, filter_sam=None, anno_file=None):
    """
    从同源基因stat结果文件挑选核心基因和特有基因
    挑选核心基因时select_sam为所有样本
    挑选单样本特有基因时，select_sam为该样本，filter_sam为其余所有样本
    挑选自定义特有基因时，select_sam的样本gene number非0，filter_sam样本gene number为0，其余样本不判断
    :param statfile: 同源基因stat结果表
    :param select_sam: 选择的样本，以逗号分隔
    :param filter_sam: 过滤的样本，以逗号分隔
    :param anno_file: 需要合并的注释表，需要有SampleGene信息
    :return:
    """
    if not os.path.exists(statfile):
        raise Exception("statfile not exists!")
    stat_table = pd.read_table(statfile, sep="\t", header=0)
    all_samples = stat_table.columns.tolist()
    if not "#ClusterID" or not "Sample_exists" in all_samples:
        raise Exception("file head has not #ClusterID or Sample_exists!")
    all_samples.remove("#ClusterID")
    all_samples.remove("Sample_exists")
    select_sams = []
    filter_sams = []
    if select_sam:
        select_sams = select_sam.split(",")
    if filter_sam:
        filter_sams = filter_sam.split(",")
    if set(all_samples) == set(select_sams):
        sam_str = ";".join(all_samples)
        select = stat_table[stat_table["Sample_exists"] == sam_str]
    elif len(select_sams) == 1 and len(filter_sams) == (len(all_samples) - 1):
        select = stat_table[stat_table["Sample_exists"] == select_sam]
    else:
        select = stat_table
        if select_sam:
            for eachs in select_sams:
                select = select[stat_table[eachs] != "-"]
        if filter_sams:
            for eachf in filter_sams:
                select = select[stat_table[eachf] == "-"]
    genes_table = contact_genes(select, select_sam, all_samples)
    if anno_file:
        result = merge_anno_file(anno_file, genes_table)
    else:
        result = genes_table.copy
    result.to_csv(outfile,sep="\t",index=False)

def contact_genes(select_table, sample, allsamples):
    if sample:
        samples = sample.split(",")
    if sample and len(samples) == 1:
        gene_table = select_table[[sample, "#ClusterID"]]
        gene_table.columns = ["Genes", "ClusterID"]
    else:
        #gene_table = pd.DataFrame(columns=("Genes", "ClusterID"))
        table_list = []
        for sam in allsamples:
            each_table = select_table[[sam, "#ClusterID"]]
            each_table.columns = ["Genes", "ClusterID"]
            table_list.append(each_table)
        gene_table = pd.concat(table_list, axis=0)
    gene_table = gene_table.drop("Genes", axis=1).join(
        gene_table["Genes"].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename(
            "Genes")).reset_index()
    tmp_sam_gene = gene_table["Genes"].str.split("|", expand=True)
    tmp_sam_gene.columns = ["Sample Name", "Gene ID"]
    tmp_cluster = gene_table["ClusterID"].str.split("(", expand=True)
    if len(tmp_cluster.columns) == 2:
        tmp_cluster.columns = ["ClusterID", "Num"]
    else:
        tmp_cluster.columns = ["ClusterID"]
    final_table = pd.concat([tmp_sam_gene, tmp_cluster["ClusterID"]], axis=1)
    final_table.drop_duplicates(['Sample Name','Gene ID'], inplace=True)  #zouguanqing 20190718
    return final_table

def merge_anno_file(anno_file, genes_table):
    anno_table = pd.read_table(anno_file, sep="\t", header=0)
    del anno_table["Gene ID"]
    if not "SampleGene" in anno_table.columns:
        raise Exception("anno_file must has column SampleGene!")
    genes_table["SampleGene"] = genes_table["Sample Name"].str.cat(genes_table["Gene ID"], sep="|")
    result = pd.merge(genes_table, anno_table, how='left', on="SampleGene")
    result = result.fillna("-")
    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, metavar="orthmcl stat file", required=True, help="orthmcl stat file")
    parser.add_argument('-anno', type=str, metavar="anno overview", required=True, help="overview anno file")
    parser.add_argument('-o', type=str, metavar="Output file", required=True, help="Output file name")
    parser.add_argument('-s', type=str, metavar="select sample", help="Input select sample names with ,")
    parser.add_argument('-f', type=str, metavar="filter sample", help="Input filter sample names with ,")
    ### select sample的样本gene number非0，filter sample样本gene number为0，其余样本不判断
    args = parser.parse_args()
    stat = args.i
    annofile = args.anno
    outfile = args.o
    select_samples = None
    filter_samples = None
    if not args.s and not args.f:
        raise Exception("must input select sample or filter sample")
    if args.s:
        select_samples = args.s
    if args.f:
        filter_samples = args.f
    select_genes(stat, outfile, select_sam=select_samples, filter_sam=filter_samples,anno_file=annofile)
