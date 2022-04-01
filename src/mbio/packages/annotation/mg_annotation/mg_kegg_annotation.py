# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang' @20190521

import sys
import pandas as pd
import argparse
from biocluster.config import Config


class MgKeggAnnotation(object):
    """
    宏基因的kegg注释线下注释方法
    """
    def __init__(self):
        self.software_dir = Config().SOFTWARE_DIR
        # self.kegg_gene = self.software_dir + '/database/KEGG/metag_database/kegg_gene_v1.xls'
        # self.kegg_ko = self.software_dir + '/database/KEGG/metag_database/kegg_ko_v1.xls'
        self.kegg_gene = self.software_dir + '/database/KEGG/metag_database/kegg_v94.2/kegg_gene_v94.2.xls'## fix by qingchen.zhang
        self.kegg_ko = self.software_dir + '/database/KEGG/metag_database/kegg_v94.2/kegg_ko_v94.2.xls'## fix by qingchen.zhang

    def run_kegg_anno(self, kegg_table, out_dir):
        """
        线下宏基因组注释流程
        :param kegg_table:kegg比对结果表m6格式；
        :param out_dir:kegg生成结果文件夹
        :return:
        """
        kegg_anno = out_dir + '/gene_kegg_anno.xls'
        kegg_enzyme = out_dir + '/kegg_enzyme_list.xls'
        kegg_module = out_dir + '/kegg_module_list.xls'
        kegg_pathway = out_dir + '/kegg_pathway_list.xls'

        align_table = pd.read_table(kegg_table, sep='\t', header=0)
        gene_table = pd.read_table(self.kegg_gene, sep='\t', header=0)
        kegg_ko = pd.read_table(self.kegg_ko, sep='\t', header=0)
        #kegg_table['Query-Name'] = (x.split('_1', 1) for x in kegg_table['Query-Name'])
        align_table['Query-Name'].replace("(_1$)", "", regex=True, inplace=True)
        tmp_gene_table = pd.merge(align_table, gene_table, left_on="Hit-Name", right_on='gene_id', how='inner')
        tmp_gene_table["web_site"] = ["http://www.genome.jp/dbget-bin/www_bget?ko:%s" % i for i in tmp_gene_table["koid"]]
        tmp_gene_table.drop_duplicates("Query-Name", "first", inplace=True)
        ko_table = pd.merge(tmp_gene_table, kegg_ko, left_on="koid", right_on='ko_id', how='inner')
        tmp_pathway_table = ko_table[['Query-Name', 'pathway_id', 'level3']]
        tmp_pathway_table = tmp_pathway_table[tmp_pathway_table['pathway_id'].astype('str') != "-"]
        tmp_enzyme_table = ko_table[['Query-Name', 'enzyme_id', 'enzyme_category']]
        tmp_enzyme_table = tmp_enzyme_table[tmp_enzyme_table['enzyme_id'].astype('str') != "-"]
        tmp_module_table = ko_table[['Query-Name', 'module_id', 'module_category']]
        tmp_module_table = tmp_module_table[tmp_module_table['module_id'].astype('str') != "-"]
        tmp_pathway_table.to_csv(kegg_pathway, sep='\t', index=False, header=0)
        tmp_enzyme_table.to_csv(kegg_enzyme, sep='\t', index=False, header=0)
        tmp_module_table.to_csv(kegg_module, sep='\t', index=False, header=0)
        tmp_anno_table = ko_table[['Query-Name', 'gene_id', 'koid', "ko_desc", 'pathway_id', 'enzyme_id', 'module_id', 'web_site','Identity-%','HSP-Len']]
        tmp_anno_table.columns = ["#Query", "Gene", "KO", "Definition", "Pathway", "Enzyme", "Module", "Hyperlink", "Identity(%)", "Align_len"]
        tmp_anno_table.to_csv(kegg_anno, sep='\t', index=False, header=1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[blast_table]',required=True,help='Input blast table')
    parser.add_argument('-o', metavar='[anno_dir]',required=True,help='output dir name')
    args = parser.parse_args()
    kegg_table = args.i
    out_dir = args.o
    meta_kegg_anno = MgKeggAnnotation()
    meta_kegg_anno.run_kegg_anno(kegg_table, out_dir)
