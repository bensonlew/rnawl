#-*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang' @20190521

import sys
import pandas as pd
import argparse
from biocluster.config import Config


class MgCogAnnotation(object):
    """
    meta宏基因cog数据详细注释信息
    """
    def __init__(self):
        self.software_dir = Config().SOFTWARE_DIR
        self.eggNOG_ID = self.software_dir + '/database/COG/eggNOG4_seqID.xls'
        self.eggNOG = self.software_dir + '/database/COG/eggNOG4.xls'

    def run_cog_anno(self, align_table, anno_table):
        """
        根据比对结果进行线下注释结果
        :param align_table: 比对输入文件 m6格式
        :param anno_table: 注释结果文件
        :return:
        """
        egg_id_table = pd.read_table(self.eggNOG_ID,sep='\t', header=1)
        egg_id_table.columns = ["nog_seq_id", "nog"]
        egg_table = pd.read_table(self.eggNOG,sep='\t', header=0)
        align_table = pd.read_table(align_table, sep='\t', header=0)
        align_table['Query-Name'].replace("(_1$)", "", regex=True, inplace=True)
        new_align_table = align_table[['Query-Name', 'Hit-Name', 'HSP-Len', 'Identity-%']]
        nog_table = pd.merge(new_align_table, egg_id_table, left_on="Hit-Name", right_on='nog_seq_id', how='inner')
        anno_table_align = pd.merge(nog_table, egg_table, on="nog", how='inner')
        #Query  NOG     NOG_description Function        Fun_description Category        Identity(%)     Align_len
        #Query-Name      Hit-Name        HSP-Len Identity-%      nog_seq_id      nog     function        category        nog_des function_des
        anno_table_align_result = anno_table_align[["Query-Name", "nog", "nog_des", "function", "function_des", "category", "Identity-%", "HSP-Len"]]
        anno_table_align_result.columns = ["#Query", "NOG", "NOG_description", "Function", "Fun_description", "Category", "Identity(%)", "Align_len"]
        anno_table_align_result.drop_duplicates("#Query", "first", inplace=True)
        anno_table_align_result.to_csv(anno_table, sep='\t', index=False, header=1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',metavar='[xml_table]',required=True,help='Input xml table')
    parser.add_argument('-o',metavar='[query_detail_table]',required=True,help='output file name')
    args = parser.parse_args()
    align_table = args.i
    anno_table = args.o
    meta_cog_anno = MgCogAnnotation()
    meta_cog_anno.run_cog_anno(align_table, anno_table)
