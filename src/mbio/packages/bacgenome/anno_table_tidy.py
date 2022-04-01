# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'

import os
import re
import pandas as pd


def anno_table_tidy(anno_table, database_type, gene_table, table_out_path):
    """

    :param anno_table: 输入注释的table文件路径
    :param database_type: 数据库类型，如nr的注释结果需特殊处理
    :param gene_table: gene预测文件，是gff格式
    :param table_out_path: 处理后注释文件
    """
    anno_table = pd.read_table(anno_table, sep='\t', header=0)
    anno_table = anno_table.ix[:, ["Query-Name", "Q-Len", "Q-Begin", "Q-End", "HSP-Len", "Hit-Description", "Hit-Name", "Hit-Len", "Hsp-Begin", "Hsp-End", "Identity-%", "E-Value", "Score"]]
    anno_table.rename(columns={'Query-Name': 'Gene ID', 'Hit-Name': 'Hit', 'Q-Len': 'Gene Len', 'Q-Begin': 'Gene Start',
                               'Q-End': 'Gene End', 'HSP-Len': 'Hit Len', 'Hsp-Begin': 'Hit Start',
                               'Hsp-End': 'Hit End', 'Identity-%': 'Identity', 'E-Value': 'Evalue', "Hit-Len":"Ref Len"}, inplace=True)
    gene_table = pd.read_table(gene_table, sep='\t', header=0)
    gene = gene_table.ix[:, ["Gene ID", "Sequence id"]]
    gene["Location"] = gene["Sequence id"].str.split("_", expand=True)[0]
    gene = gene.ix[:, ["Gene ID", "Location"]]
    if database_type == "nr":
        anno_table["Hit-Description"] = list(map(nr_anno_des, anno_table["Hit-Description"]))
        table_out = gene.merge(anno_table, on='Gene ID', how='left')
        table_out = table_out.fillna("-")
        #table_out["Hit-Description"] = [c.replace("^-$", "hypothetical protein") for c in table_out["Hit-Description"]]
        table_out["Hit-Description"] = [re.subn("^-$", "hypothetical protein", c)[0] for c in table_out["Hit-Description"]]
        table_out = table_out[~ table_out['Hit'].isin(['-'])]
    else:
        table_out = gene.merge(anno_table, on='Gene ID', how='left')
        table_out = table_out.fillna("-")
    table_out.index = table_out["Gene ID"]
    table_out =table_out.ix[:,["Location", "Hit","Ref Len", "Hit-Description", "Gene Len", "Gene Start", "Gene End", "Hit Len", "Hit Start",	"Hit End", "Identity", "Evalue", "Score"]]
    table_out.to_csv(table_out_path, sep="\t")
    return table_out

def nr_anno_des(des):
    """
    描述信息处理：
    1. nr的注释信息要去掉“MULTISPECIES: ”以及最后中括号和里面的物种信息，如MULTISPECIES: oligoribonuclease oligoribonuclease [Vibrio sp. HENC-01] 最后变成“oligoribonuclease oligoribonuclease”。
    2. 如果注释信息里面包括“partial”，如hypothetical protein, partial possible integrase, partial [Lactobacillus fermentum]等，都要去掉“,partial"和物种信息。
    """
    des = re.subn("^ MULTISPECIES: ", " ", des)[0]
    des = re.subn(" \[.*\]", "", des)[0]
    des = re.subn(", partial.*", "", des)[0]
    des = re.subn("[hH]ypothetical protein.*", "hypothetical protein", des)[0]
    return des

def get_nr_des(anno_table, nr_table, table_out_path,database=None):
    anno_table = pd.read_table(anno_table, sep='\t', header=0)
    anno_table = anno_table.rename(columns={'Seq_id':'#Query'})
    nr_table = pd.read_table(nr_table, sep='\t', header=0)
    nr = nr_table.ix[:, ["Gene ID", "Location", "Hit-Description"]]
    nr.columns = ["Gene ID", "Location", "Gene Description"]
    nr = nr.drop_duplicates(subset=["Gene ID", "Location"], keep='first')
    table_out = anno_table.merge(nr, left_on='#Query', right_on = 'Gene ID', how='left')
    table_out.dropna(subset=['Gene ID'], inplace=True)
    table_out.index = table_out["Gene ID"]
    del table_out["Gene ID"]
    del table_out["#Query"]
    if database == "cazy":
        table_out["Coverage(%)"] = table_out["Coverd_fraction"]*100
        table_out = table_out["Gene ID", "Location", "Family", "Family Description", "Class","Class Description","Identity(%)","Coverage(%)","Evalue","Score"]

    table_out.to_csv(table_out_path, sep="\t")

if __name__ == '__main__':  # for test
    #anno_table = "/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/GH01/RefalignDnabac/YC_YT1_predict_vs_CP010953.1.faa.xls"
    anno_table = "/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/GH01/NR/blast.table"
    #anno_table = "/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/GH01/Swissprot/swissprot.xls"
    #anno_table = "/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/GH01/Cazy/gene_cazy_parse_anno.xls"
    database_type = "nr" #"ref"
    gene_table= "/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/GH01/Predict/YC_YT1_predict.gff"
    #table_out_path= "/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/GH01/RefalignDnabac/anno_ref.xls"
    #table_out_path= "/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/GH01/NR/anno_nr.xls"
    #table_out_path= "/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/GH01/Swissprot/anno_swissprot.xls"
    #table_out_path= "/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/GH01/Cazy/anno_cazy.xls"
    #anno_table_tidy(anno_table, database_type, gene_table, table_out_path)
    nr = "/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/GH01/NR/anno_nr.xls"
    cog = "/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/GH01/COG/gene_cog_anno.xls"
    kegg = "/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/GH01/KEGG/gene_kegg_anno.xls"
    go = "/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/GH01/GO/query_gos.list"
    cogpath= "/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/GH01/COG/cog_anno.xls"
    cazy = "/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/GH01/Cazy/gene_cazy_parse_anno.xls"
    keggpath = "/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/GH01/KEGG/kegg_anno.xls"
    gopath = "/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/GH01/GO/anno_go.xls"
    cazypath = "/mnt/ilustre/users/sanger-dev/sg-users/zhujuan/Genomic/GH01/Cazy/anno_cazy.xls"
    #get_nr_des(cazy,nr,cazypath)

