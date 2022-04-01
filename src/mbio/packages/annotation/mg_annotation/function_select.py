# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

import sys
import argparse
import os
import pandas as pd
from mbio.packages.statistical.large_data import table_parser


class mg_function_select(object):
    def __init__(self):
        self.levels_keep = {}  # 保留的层级
        self.levels_del = {}  # 去除的层级
        self.database = ''
        self.lowest_level = {
            "nr": "Species",
            "cog": "NOG",
            "kegg": "Gene",
            "ardb": "ARG",
            "card": "ARO",
            "vfdb": "VFs",
            "cazy": "Family",
            "probio": "Probiotic_name",
            "go": "GO Term (Lev4)",
            "phi": "protein",
            "mvirdb": "Virulence Factor ID",
            "qs": "QS_id",
            "pfam": "Pfam ID",
            "tcdb": "TCDB ID",
            "p450": "Sid"
        }

    def read_select_level(self, level_name_file):
        with open(level_name_file, "r") as f:
            for line in f:
                line = line.strip().split("\t")
                level = line[1]
                level_name = line[2]
                if int(line[0]) == 0:
                    self.levels_del.setdefault(level, []).append(level_name)
                elif int(line[0]) == 1:
                    self.levels_keep.setdefault(level, []).append(level_name)
                else:
                    raise Exception("select_level file error!")

    def from_gene_select_anno(self, gene_list, anno_file):
        reader1 = pd.read_table(gene_list, sep='\t', header=0, iterator=True)
        gene_list = table_parser(reader1, chunk_size=100000, ignore_index=True)
        genes = gene_list['GeneID']
        select_anno = anno_file[anno_file['#Query'].isin(genes)]
        return select_anno

    def creat_new_anno(self, annofile, level_name_file, database, iden, align_len, outfile, lowest_level =None,
                       gene_list=None):
        self.database = database
        reader = pd.read_table(annofile, sep='\t', header=0, iterator=True)
        anno_file = table_parser(reader, chunk_size=100000, ignore_index=True).drop_duplicates('#Query')
        if gene_list != None:
            anno_file = self.from_gene_select_anno(gene_list, anno_file)
            '''
            if len(anno_file) < 1:
                raise Exception("该基因集注释表结果为空！")
            '''
        if lowest_level != None:
            lowest_level_name = self.lowest_level[database]
            reader2 = pd.read_table(lowest_level, sep='\t', header=0, iterator=True)
            sel_lowest = table_parser(reader2, chunk_size=100000, ignore_index=True)
            sel_lowest = pd.DataFrame(sel_lowest)
            lowest_names = sel_lowest['#name']
            if database == "nr":
                lowest_names = lowest_names.apply(lambda x: x.split(";")[-1])
            anno_file = anno_file[anno_file[lowest_level_name].isin(lowest_names)]
            '''
            if len(anno_file) < 1:
                raise Exception("该基因集下丰度筛选结果为空！")
            '''
        print anno_file.head()
        if iden != 0:
            if 'Probability(%)' in anno_file.columns:
                anno_file = anno_file[(anno_file['Probability(%)'] >= float(iden) * 100)]
            else:
                anno_file = anno_file[(anno_file['Identity(%)'] >= float(iden) * 100)]
        if align_len != 0:
            anno_file = anno_file[anno_file['Align_len'] >= int(align_len)]
        '''
        anno_file = anno_file[
            (anno_file['Identity(%)'] >= float(iden) * 100) & (anno_file['Align_len'] >= int(align_len))]
        if len(anno_file) < 1:
            raise Exception("identity或length筛选结果为空！")
        '''
        if level_name_file == "all":
            anno_file.to_csv(outfile, sep="\t", index=False)
        else:
            self.read_select_level(level_name_file)
            choose_tables = []
            if len(self.levels_keep.keys()) != 0:
                for k_level in self.levels_keep.keys():
                    print k_level
                    levelnames = self.levels_keep[k_level]

                    tmp_choose = self.select(anno_file, k_level, levelnames, database, 0)
                    choose_tables.append(tmp_choose)
                choose = pd.concat(choose_tables)
            else:
                choose = anno_file
            if len(self.levels_del.keys()) != 0:
                for d_level in self.levels_del.keys():
                    levelnames = self.levels_del[d_level]
                    choose = self.select(choose, d_level, levelnames, database, 1)
            '''
            if len(choose) < 1:
                raise Exception("功能筛选结果为空！")
            '''
            choose = self.revert(choose)
            choose.to_csv(outfile, sep="\t", index=False)

    def select(self, table, level, choose_names, database, type):
        '''
         type 0 保留，1 去除
         不仅筛选基因行并修改该行内容
         return select_table
        '''
        new_table = table.copy()
        if level == 'Class' and self.database == 'cazy':
            if 'Class_description' in new_table.columns:
                print('Class changed to Class_description')
                level = 'Class_description'
        level_choose = new_table[level].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename(level)
        level_choose = pd.DataFrame(level_choose)
        level_choose.columns = [level]
        interlock_level = ["Category", "Function","Pathway", "Level1", "Level2", "Level3", 'Class']
        if database in ["cog","kegg", 'vfdb', 'card'] and level in interlock_level:
            level_choose = self.split_level(database, level, level_choose, table)
        if type == 0:  # 0 保留，1 去除
            new = level_choose[level_choose[level].isin(choose_names)]
        else:
            new = level_choose[~level_choose[level].isin(choose_names)]
        result = self.revert(new)
        new_table[result.columns] = result
        new_table= new_table.dropna()
        return new_table

    def split_level(self, database, level, level_choose, table):
        '''
         database :注释数据库，cog和kegg时用
         level: 筛选的层级
         level_choose: 分号分隔并用stack堆叠的筛选层级table
         table:原始table
         return：筛选cog或kegg的某些层级时，其他层级联合改动的level堆叠表
        '''
        rep_level = []
        if database == "cog":
            rep_level = ["Category", "Fun_description", "Function"]
        if database == "kegg":
            rep_level = ["Pathway", "Level1", "Level2", "Level3"]
        if database == "vfdb":
            rep_level = ["Level1", "Level2"]
        if database == "card":
            rep_level = ["Class", "Class_description"]
        rep_level.remove(level)
        print rep_level
        for each in rep_level:
            level_table = table[each].str.split(';', expand=True).stack().reset_index(level=1, drop=True)
            level_table = pd.DataFrame(level_table)
            level_table.columns = [each]
            level_choose = pd.concat([level_choose, level_table], axis=1)
        return level_choose

    def revert(self, level_table):
        '''
         level_table: 筛选过后的table
         return：相同index的行以分号重新连接的table
        '''
        levels = ",".join(level_table.columns)
        level_list = levels.split(",")
        level_table = pd.DataFrame(level_table).astype(str)
        for i in range(0, len(level_list)):
            if i == 0:
                result = level_table.groupby(level_table.index)[level_list[i]].apply(';'.join)
                revert_table = pd.DataFrame(result)
                revert_table.columns = [level_list[i]]
            else:
                tmp_result = level_table.groupby(level_table.index)[level_list[i]].apply(';'.join)
                tmp_result = pd.DataFrame(tmp_result)
                tmp_result.columns = [level_list[i]]
                revert_table = revert_table.join([tmp_result])
        return revert_table


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='[gene_*_anno.xls]', required=True, help='Input gene annotation table')
    parser.add_argument('-g', metavar='[geneset file]', help='Input gene set table')
    parser.add_argument('-s', metavar='[select lowest file]', help='Input select lowest function table')
    parser.add_argument('-o', metavar='[output file]', required=True, help='output file name')
    parser.add_argument('-l', metavar='[select level file]', help='input select level file', default="all")
    parser.add_argument('-iden', metavar='identity', help='input  identity cutoff', default=0.0)
    parser.add_argument('-len', metavar='align_length', help='input align_length cutoff', default=0)
    parser.add_argument('-database', metavar='[database]', required=True, help='input database name')
    args = parser.parse_args()
    annofile = args.i
    level_name_file = args.l
    iden = args.iden
    align_len = args.len
    database = args.database
    outfile = args.o
    function_select = mg_function_select()
    if args.s:
        lowest_level = args.s
    else:
        lowest_level = None
    if args.g:
        gene_list = args.g
    else:
        gene_list = None
    function_select.creat_new_anno(annofile, level_name_file, database, iden, align_len, outfile,
                                   lowest_level=lowest_level, gene_list=gene_list)

