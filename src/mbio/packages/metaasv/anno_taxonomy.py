# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang' @20200511

import pandas as pd
import argparse
from biocluster.config import Config


class AnnoTaxon(object):
    """
    注释
    """
    def __init__(self):
        self.software_dir = Config().SOFTWARE_DIR

    def run_taxon_anno(self, blast_table, taxonomy, out_file):
        """
        blast方法注释流程
        根据线下注释的tax得到结果
        :param blast_table:nr或nt比对结果表m6格式的比对结果表；
        :param taxonomy:参考注释信息
        :param out_dir:比对结果
        :return:
        """
        taxon_table = pd.read_table(taxonomy, sep='\t', header=None).astype(str)
        align_table = pd.read_table(blast_table, sep='\t', header=None).astype(str)
        tmp_gene_table = pd.merge(align_table, taxon_table, left_on=1, right_on=0, how='inner')
        tmp_gene_table.columns = ["Subject_ID1","ASV_ID", "Subject_ID", "Identity", "Alignment_length", "Mismatches", "Gap", "Query_start", "Query_end", "Subject_start", "Subject_end", "E-value", "Score","Subject_ID2", "Taxonomy"]
        tmp_anno_table = tmp_gene_table[["ASV_ID", "Taxonomy", "Identity"]]
        tmp_anno_table["Identity"] = tmp_anno_table["Identity"].apply(lambda x:float(x)/100)
        tmp_anno_table.to_csv(out_file, sep='\t', index=False, header=0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[xml_table]',required=True,help='Input xml table')
    parser.add_argument('-a', metavar='[database]', required=True, help='type of taxon db')
    parser.add_argument('-o', metavar='[anno_file]',required=True,help='output file name')
    args = parser.parse_args()
    kegg_table = args.i
    database = args.a
    out_file = args.o
    taxon_anno = AnnoTaxon()
    taxon_anno.run_taxon_anno(kegg_table, database, out_file)