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
        self.hmdb_ref = Config().SOFTWARE_DIR + "/database/HMDB/hmdb_level.xls"

    def creat_anno(self, overviewfile, outDir, metabset=None):
        """
        HMDB注释和统计主函数
        :param table: 代谢物注释详情表
        :param outDir: 输出目录
        :param metabset: 代谢集文件, 没有head
        :return: 代谢物hmdb注释表文件
        """
        if not os.path.exists(outDir):
            os.mkdir(outDir)
        print self.hmdb_ref
        ref_hmdb = pd.read_table(self.hmdb_ref, sep="\t", header=0)
        overview_table = pd.read_table(overviewfile, sep="\t", header=0)
        name_order = overview_table.columns
        # overview_table1 = overview_table.drop("hmdb_id", axis=1).join(
        #     overview_table["hmdb_id"].str.split(";", expand=True).stack().reset_index(
        #         level=1, drop=True).rename("hmdb_id"))
        # overview_table1 = overview_table1[name_order]
        def deal_hmdb(x):
            spx = x.split(';')
            for e in spx:
                if re.match('^HMDB\d+$', e):
                    return e
            return x

        metab_hmmdb = overview_table[['metab_id','hmdb_id']]
        metab_hmmdb['hmdb_id'] = metab_hmmdb['hmdb_id'].apply(lambda x: deal_hmdb(x))
        tmp_merge = pd.merge(metab_hmmdb,ref_hmdb,how='left',on='hmdb_id')
        tmp_merge.drop('hmdb_id',axis=1)
        anno_table = pd.merge(overview_table,tmp_merge, how='left', on='metab_id')

        #anno_table = pd.merge(overview_table, ref_hmdb, how='left', on="hmdb_id")
        anno_table = anno_table.fillna("-")
        outoverview = os.path.join(outDir, "anno.xls")
        tmp_table = anno_table[["kingdom","super_class","class","sub_class"]]
        tmp_table.index = anno_table["metab_id"]
        overview_table.index = overview_table["metab_id"]
        outoverview_table = pd.merge(overview_table, tmp_table, how='left', left_index=True, right_index=True, sort=False)
        outoverview_table = outoverview_table.drop_duplicates("metab_id")
        outoverview_table.to_csv(outoverview, sep="\t",index=False)
        anno_table = anno_table[["metab", "metab_id", "kingdom", "super_class", "class", "sub_class"]]
        anno_table.columns = ["Metabolite", "Metab_id", "Kingdom", "Superclass", "Class", "Subclass"]
        level_file = os.path.join(outDir, "HmdbLevel_Origin.xls")
        anno_table = anno_table[anno_table["Class"] != "-"]
        anno_table = anno_table.drop_duplicates("Metab_id")
        anno_table.to_csv(level_file, sep="\t", index=False)
        if metabset:
            metabset_table = pd.read_table(metabset, sep="\t", header=0)
            anno_table = pd.merge(anno_table, metabset_table, how='inner', left_on="Metab_id",right_on="metab_id")
            level_file2 = os.path.join(outDir, "HmdbLevel.xls")
            anno_table = anno_table[["Metabolite", "Metab_id", "Kingdom", "Superclass", "Class", "Subclass"]]
            anno_table.to_csv(level_file2, sep="\t", index=False)
            return level_file2
        return level_file

    def select_anno(self, hmdb_anno, select_file, outDir, header="F"):
        """
        :param table: HMDB注释
        :param outDir: 输出目录
        :param metabset: 代谢集文件, 没有head
        :return: 代谢物hmdb注释表文件
        """
        if not os.path.exists(outDir):
            os.mkdir(outDir)
        if header == "T":
            select_table = pd.read_table(select_file, sep="\t", header=0)
        else:
            select_table = pd.read_table(metabset, sep="\t", header=None)
            select_table.columns = ["metab_id"]
        hmdb_anno = pd.read_table(hmdb_anno, sep="\t", header=0)
        anno_table = pd.merge(hmdb_anno, select_table, how='inner', left_on="Metab_id", right_on="metab_id")
        del anno_table['metab_id']
        level_file = os.path.join(outDir, "HmdbLevel.xls")
        anno_table.to_csv(level_file, sep="\t", index=False)
        return level_file

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
        subclass_stat = dict(list(anno_table["Metab_id"].groupby([anno_table[levelname]])))
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

    def run_main(self, overview, outDir, metabset=None, origin="F", header="F"):
        if origin == "T":
            '''
            工作流使用，生成原始的所有anno_overview生成的hmdb_level_origin.xls
            以及跟阴阳离子和mix表相关的hmdb_level.xls
            '''
            print "anno"
            level_file = self.creat_anno(overview, outDir, metabset=metabset)
            self.creat_level(level_file, outDir)
        else:
            '''
            交互分析使用
            '''
            print "select"
            level_file = self.select_anno(overview, metabset, outDir, header=header)
            self.creat_level(level_file, outDir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, metavar="overview anno", required=True, help="overview anno file")
    parser.add_argument('-o', type=str, metavar="Output directory", required=True, help="Output directory name")
    parser.add_argument('-s', type=str, metavar="metabset", help="Input metabset file")
    parser.add_argument('--origin', action='store_true', help="is for origin workflow table")
    parser.add_argument('--header', action='store_true', help="is selcet table with header")
    parser.add_argument('-mongo', help="mongo verison", type=int, default=None)
    args = parser.parse_args()
    if args.mongo:
        Config().DBVersion = args.mongo
    overview = args.i
    outDir = args.o
    origin = "F"
    header = "F"
    metabset = None
    if args.s:
        metabset = args.s
    if args.origin:
        origin = "T"
    if args.header:
        header = "T"
    run = AnnoHmdb()
    run.run_main(overview, outDir, metabset=metabset, origin=origin, header=header)
