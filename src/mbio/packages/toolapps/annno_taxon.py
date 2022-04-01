# -*- coding: utf-8 -*-
# __author__ = 'hao.gao' @20200119

import pandas as pd
import argparse
from biocluster.config import Config


class AnnoTaxon(object):
    """
    NR或NT的物种注释信息
    """
    def __init__(self):
        self.software_dir = Config().SOFTWARE_DIR
        self.nr_taxon = self.software_dir + '/database/toolapps/NR_2018.06.06/taxon/NR_tax.xls'
        self.nt_taxon = self.software_dir + '/database/toolapps/NT_2018.07.26/taxon/NT_tax.xls'
        self.other = self.software_dir + '/../workspace/20200116/Single_tsg_36852/Blast/ReadsAlignDatabase/SampleFq/FqReads/ReadsAlign/Blast/output/bb'

    def run_taxon_anno(self, blast_table, database, out_file):
        """
        线下细菌因组注释流程
        :param blast_table:nr或nt比对结果表m6格式；
        :param database:比对数据库类型
        :param out_dir:比对数据库类型
        :return:
        """
        taxon_table = ''
        if database in ["NR"]:
            taxon_table = pd.read_csv(self.nr_taxon, sep='\t', header=0, chunksize=100000)
        elif database in ["NT"]:
            taxon_table = pd.read_csv(self.nt_taxon, sep='\t', header=0, chunksize=100000)
        else:
            taxon_table = pd.read_csv(self.other, sep='\t', header=0)
        align_table = pd.read_table(blast_table, sep='\t', header=None)
        align_table.columns =["q_id","assession","c","d","e","f","g","h","i","j","k","l"]
        list1= []
        for chunk in taxon_table:
            tmp_gene_table = pd.merge(align_table, chunk, left_on="assession", right_on='accession.version_number', how='inner')
            tmp_anno_table = tmp_gene_table[["q_id", 'taxonomy']]
            list1.append(tmp_anno_table)
        df = pd.concat(list1,ignore_index = True)
        df.columns = ["#Sequence", "Taxonomy"]
        df.to_csv(out_file, sep='\t', index=False, header=1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[xml_table]',required=True,help='Input xml table')
    parser.add_argument('-a', metavar='[database]', required=True, help='type of database')
    parser.add_argument('-o', metavar='[anno_file]',required=True,help='output file name')
    args = parser.parse_args()
    kegg_table = args.i
    database = args.a
    out_file = args.o
    taxon_anno = AnnoTaxon()
    taxon_anno.run_taxon_anno(kegg_table, database, out_file)
