# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'


import argparse
from biocluster.config import Config
import pandas as pd


class meta_vfdb_anno(object):
    """
    meta宏基因vfdb数据详细注释信息
    注释数据库版本vfdb_20200703
    改为线下注释
    """
    def __init__(self):
        self.software_dir = Config().SOFTWARE_DIR
        self.vfdb_database = self.software_dir + '/database/VFDB/vfdb_v20200703_database.xls'

    def run_vfdb(self, align_table, anno_table):
        """
        根据比对结果表和数据库得到注释信息
        :param align_table: 输入比对结果表
        :param anno_table: 输出注释结果表
        :return:
        """
        database_table = pd.read_table(self.vfdb_database,sep='\t', header=1)
        ##vfg_id  gi_number       gene    gene_des        vfs_function    vfs     origin  level1  level2
        align_table = pd.read_table(align_table, sep='\t', header=0)
        ##Score E-Value HSP-Len Identity-%  Similarity-%  Query-Name  Q-Len  Q-Begin Q-End  Q-Frame Hit-Name  Hit-Len Hsp-Begin Hsp-End Hsp-Frame Hit-Description
        align_table['Query-Name'].replace("(_1$)", "", regex=True, inplace=True)
        new_align_table = align_table[['Query-Name', 'Hit-Name', 'HSP-Len', 'Identity-%']]
        new_align_table['Hit-Name'] = new_align_table['Hit-Name'].str.split("(")[0]
        vfdb_table = pd.merge(new_align_table, database_table, left_on="Hit-Name", right_on='vfg_id', how='inner')
        anno_table_align = vfdb_table[['Query-Name',"gi_number", "gene", "gene_des", "vfs", "vfs_function", "origin", "level1", "level2", 'Identity-%','HSP-Len']]
        anno_table_align.columns = ['Query-Name',"Gi_number", "VfdbGene", "Gene_description", "VFs", "VF_Function", "Species", "Level1", "Level2", 'Identity(%)','Align_len']
        anno_table_align.to_csv(anno_table, sep='\t', index=False, header=1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',metavar='[xml_table]',required=True,help='Input xml table')
    parser.add_argument('-o',metavar='[query_detail_table]',required=True,help='output file name')
    args = parser.parse_args()
    align_table = args.i
    anno_table = args.o
    meta_vfdb_anno = meta_vfdb_anno()
    meta_vfdb_anno.run_vfdb(align_table, anno_table)

