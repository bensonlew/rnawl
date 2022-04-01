# -*- coding: utf-8 -*-
import pandas as pd
from biocluster.api.database.base import Base
from biocluster.config import Config
import os
import argparse
import re


class AnnoHmdb(Base):
    def __init__(self):
        super(AnnoHmdb, self).__init__()
        self._project_type = 'metabolome'
        self.hmdb_ref = Config().SOFTWARE_DIR + "/database/HMDB/hmdb_metabolites.detail.xls"

    def creat_anno(self, overviewfile, outDir):
        """
        HMDB注释和统计主函数
        :param table: 代谢物注释详情表
        :param outDir: 输出目录
        :param metabset: 代谢集文件, 没有head
        :return: 代谢物hmdb注释表文件
        """
        if not os.path.exists(outDir):
            os.mkdir(outDir)
        print(self.hmdb_ref)
        ref_hmdb = pd.read_table(self.hmdb_ref, sep="\t", header=0)
        overview_table = pd.read_table(overviewfile, sep="\t", header=0)


        metab_hmmdb = overview_table[['metab_id','name']]

        anno_table = pd.merge(metab_hmmdb,ref_hmdb,how='left',on='name')
        #tmp_merge.drop('hmdb_id',axis=1)

        anno_table = anno_table.fillna("-")


        anno_table = anno_table[["metab_id", "name","accession", "kingdom", "super_class", "class", "sub_class"]]
        anno_table.columns = ["metab_id","Metabolite","HMDB ID", "Kingdom", "Superclass", "Class", "Subclass"]
        detail_file = os.path.join(outDir, "Hmdb_detail.xls")
        #####anno_table = anno_table[anno_table["Class"] != "-"]
        #anno_table = anno_table.drop_duplicates("Metab_id")
        anno_table.to_csv(detail_file, sep="\t", index=False)

        return detail_file



    def creat_level(self, anno_file, outDir):
        """
        统计level层级代谢物count数
        :param anno_file: 代谢物注释详情表文件
        :param outDir: 输出目录
        :return:
        """
        mydict = {}
        if os.path.exists(anno_file):
            anno_table = pd.read_table(anno_file, sep="\t", header=0)
        else:
            raise Exception("anno_file pat error")
        if len(anno_table) >0:
            level_list = ["Superclass", "Class", "Subclass"]
            for eachlevel in level_list:
                result_table = self.creat_each_level(anno_table, eachlevel)
                outfile = os.path.join(outDir, "Hmdb" + eachlevel + ".xls")
                result_table.to_csv(outfile, sep="\t", index=False)

    def creat_each_level(self, anno_table, levelname):
        mydict = {}
        subclass_stat = dict(list(anno_table["metab_id"].groupby([anno_table[levelname]])))
        for each in subclass_stat.keys():
            metab_ids = subclass_stat[each].tolist()
            number = len(metab_ids)
            metab_ids = ";".join(metab_ids)
            mydict[each] = [number, metab_ids]
        result = pd.DataFrame(mydict).T
        result.columns = ["Number", "Metab_ids"]
        result[levelname] = result.index
        result = result[[levelname, "Number", "Metab_ids"]]
        result = result[(result[levelname] != "-")]
        result = result.sort_values(["Number"], ascending=False)
        return result


    def run_main(self, overview, outDir):
        level_file = self.creat_anno(overview, outDir)
        self.creat_level(level_file, outDir)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, metavar="metab list", required=True, help="metab list")
    parser.add_argument('-o', type=str, metavar="Output directory", required=True, help="Output directory name")

    args = parser.parse_args()

    overview = args.i
    outDir = args.o

    run = AnnoHmdb()
    run.run_main(overview, outDir)
