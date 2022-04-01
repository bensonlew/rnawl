# -*- coding: utf-8 -*-
# __author__ = 'hao.gao' @20190919

import sys
import pandas as pd
import argparse
from biocluster.config import Config


class KeggAnno(object):
    """
    细菌基因组的kegg注释线下注释方法
    """
    def __init__(self):
        self.software_dir = Config().SOFTWARE_DIR
        self.kegg_gene = self.software_dir + '/database/bac_kegg/kegg_function.xls'

    def run_kegg_anno(self, kegg_table, out_dir):
        """
        线下细菌因组注释流程
        :param kegg_table:kegg比对结果表m6格式；
        :param out_dir:kegg生成结果文件夹
        :return:
        """
        kegg_anno = out_dir + '/gene_kegg_anno.xls'
        align_table = pd.read_table(kegg_table, sep='\t', header=0)
        gene_table = pd.read_table(self.kegg_gene, sep='\t', header=0)
        tmp_gene_table = pd.merge(align_table, gene_table, left_on="Hit-Name", right_on='gene_id', how='inner')
        tmp_anno_table = tmp_gene_table[['Query-Name', 'Gene_name', 'KO', "KO_des", 'Ko_id', 'Enzyme', 'Module', 'Hyperlink', 'Level1', 'Level2', 'Level3', 'HSP-Len', 'Identity-%']]
        tmp_anno_table.columns = ["#Query", "Gene Name", "KO", "Definition", "Pathway", "Enzyme", "Module", "Hyperlink", 'Level1', 'Level2', 'Level3', "Identity(%)", "Align_len"]
        tmp_anno_table.to_csv(kegg_anno, sep='\t', index=False, header=1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[xml_table]',required=True,help='Input xml table')
    parser.add_argument('-o', metavar='[anno_dir]',required=True,help='output dir name')
    args = parser.parse_args()
    kegg_table = args.i
    out_dir = args.o
    meta_kegg_anno = KeggAnno()
    meta_kegg_anno.run_kegg_anno(kegg_table, out_dir)
