# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'


import argparse
from biocluster.config import Config
import pandas as pd


class meta_card_anno(object):
    """
    meta宏基因card数据详细注释信息
    注释数据库版本card_v3.0.9
    改为线下注释
    """
    def __init__(self):
        self.software_dir = Config().SOFTWARE_DIR
        self.card_database = self.software_dir + '/database/card/card_v3.0.9_aro_index.tsv'
        self.aro_des = self.software_dir + '/database/card/card_v3.0.9_aro.tsv'

    def run_card(self, align_table, anno_table):
        """
        根据比对结果表和数据库得到注释信息
        :param align_table: 输入比对结果表
        :param anno_table: 输出注释结果表
        :return:
        """
        with open(align_table, 'r') as infile, open(anno_table, 'wb') as outfile:
            outfile.write('#Query\tARO\tARO_name\tARO_description\tARO_category\tClass\tClass_description\tIdentity(%)\tAlign_len\n')
        database_table = pd.read_table(self.card_database,sep='\t', header=0)
        ##ARO Accession	CVTERM ID Model Sequence ID	Model ID Model Name	ARO Name Protein Accession	DNA Accession AMR Gene Family Drug Class	Resistance Mechanism
        align_table = pd.read_table(align_table, sep='\t', header=0)
        ##Score E-Value HSP-Len Identity-%  Similarity-%  Query-Name  Q-Len  Q-Begin Q-End  Q-Frame Hit-Name  Hit-Len Hsp-Begin Hsp-End Hsp-Frame Hit-Description
        align_table['Query-Name'].replace("(_1$)", "", regex=True, inplace=True)
        new_align_table = align_table[['Query-Name', 'Hit-Name', 'HSP-Len', 'Identity-%']]
        card_table = pd.merge(new_align_table, database_table, left_on="Hit-Name", right_on='ARO Accession', how='inner')
        card_columns = list(card_table.columns)
        card_columns.remove("Description")
        card_table = card_table[card_columns]
        aro_table = pd.read_table(self.aro_des,sep='\t', header=0)
        ##Accession       Name    Description     ID
        all_anno_table = pd.merge(card_table, aro_table, left_on="ARO Accession", right_on='Accession', how='left')
        anno_table_align = all_anno_table[['Query-Name',"Accession", "Name", "Description", "AMR Gene Family", "Drug Class", "Type","Resistance Mechanism", 'Identity-%','HSP-Len']]
        anno_table_align.columns = ['#Query',"ARO", "ARO_name", "ARO_description", "AMR_Gene_Family", "Drug_Class","Antibiotic_class", "Resistance_Mechanism", 'Identity(%)','Align_len']
        anno_table_align.to_csv(anno_table, sep='\t', index=False, header=1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',metavar='[xml_table]',required=True,help='Input xml table')
    parser.add_argument('-o',metavar='[query_detail_table]',required=True,help='output file name')
    args = parser.parse_args()
    align_table = args.i
    anno_table = args.o
    meta_card_anno = meta_card_anno()
    meta_card_anno.run_card(align_table, anno_table)

