#-*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang' @20191226


import pandas as pd
import argparse
import os
import subprocess
from biocluster.config import Config


class PickAnnotation(object):
    """
    根据输入otu表，得到注释信息
    """
    def __init__(self):
        self.software_dir = Config().SOFTWARE_DIR
        kegg_database = self.software_dir + '/database/meta/picrust2/KEGG/kegg_v94.2_database.xls'
        cog_database = self.software_dir + '/database/meta/picrust2/COG/cog_database.xls'
        metacyc_database = self.software_dir + '/database/meta/picrust2/Metacyc/metacyc_database.txt'
        # kegg_database = '/mnt/ilustre/users/sanger-dev/home/zhangqingchen/metagenomic/picrust/KEGG/kegg_database.xls'
        # cog_database = '/mnt/ilustre/users/sanger-dev/home/zhangqingchen/metagenomic/picrust/COG/cog_database.xls'
        # metacyc_database = '/mnt/ilustre/users/sanger-dev/home/zhangqingchen/metagenomic/picrust/Metacyc/metacyc_database.txt'
        self.cog_func_database = self.software_dir + "/database/meta/picrust2/COG/function.xls"
        self.ko_database = self.software_dir + "/database/meta/picrust2/KEGG/kegg_v94.2_level_databse.xls"
        self.ec_database = self.software_dir + "/database/meta/picrust2/KEGG/ec_level4_info.tsv"
        # self.module_database = self.software_dir + "/database/meta/picrust2/KEGG/pick_Moule.xls"
        self.module_database = self.software_dir + "/database/meta/picrust2/KEGG/pick_Moule_v94.2.xls" ## module的新库
        self.database = {
            "cog": cog_database,
            "kegg": kegg_database,
            "metacyc": metacyc_database,
        }

    def run_cog_anno(self, otu_table, out_dir, database):
        """
        根据输入的otu表得到cog的注释结果
        :param otu_table: 输入的otu丰度表
        :param out_dir: 输出文件夹
        :return:
        """
        print("开始进行cog注释分析")
        otu_path = otu_table
        otu_name = os.path.basename(otu_path)
        cog_outfile = os.path.join(out_dir, "prediction_cog.xls")
        cog_table = 'cog_abundance.xls'
        function_outfile = os.path.join(out_dir, "prediction_function.xls")
        if os.path.exists(cog_outfile):
            os.remove(cog_outfile)
        if os.path.exists(function_outfile):
            os.remove(function_outfile)
        if otu_name.endswith(".gz"):
            otu_file = self.unzip_file(otu_path, cog_table)
        else:
            otu_file = otu_table
        data = pd.read_table(otu_file, sep="\t", header=0)
        sample_list = list(data.columns)[1:]
        data.columns = ["cog_id"] + sample_list
        database_file = self.database[database]
        data_database = pd.read_table(database_file, sep="\t", header=0)
        merge_file = pd.merge(data_database, data, left_on='nog', right_on='cog_id', how='inner')
        cog_file_list = ["cog_id", "nog_des", "category"] + sample_list
        function_list = ["function"] + sample_list
        cog_file = merge_file.loc[:, cog_file_list]
        cog_file.drop_duplicates("cog_id", "first", inplace=True)
        cog_file.columns = ["COG_ID", "COG_Description", "Function"] + sample_list
        cog_file.to_csv(cog_outfile, sep='\t', index=False, header=1)
        function_file = merge_file.loc[:, function_list]
        function_file2 = function_file.drop("function", axis=1).join(function_file['function'].str.split(";", expand=True).stack().reset_index(level=1,drop=True).rename('function'))
        function_file3 = function_file2.groupby("function").sum()
        # function_file3.drop_duplicates("function", "first", inplace=True)
        print(function_file3.columns)
        function_file3 = function_file3.reset_index()
        function_file3.columns = ["function"] + sample_list
        function_database = pd.read_table(self.cog_func_database, sep="\t", header=0)
        function_merge_file = pd.merge(function_database, function_file3, left_on='function', right_on='function', how='inner')
        function_outfile_list = ["function", "function_description"] + sample_list
        function_out = function_merge_file.loc[:, function_outfile_list]
        function_out.columns = ["Function", "Function_Description"] + sample_list
        function_out.to_csv(function_outfile, sep='\t', index=False, header=1)

    def run_kegg_anno(self, otu_table, out_dir, database, enzyme):
        """
        根据输入的otu表得到cog的注释结果
        :param otu_table: 输入的KO的otu丰度表
        :param out_table: 输出结果丰度表
        :param enzyme: 输出结果丰度表
        :return:
        """
        print("开始进行kegg注释分析")
        otu_path = otu_table
        otu_name = os.path.basename(otu_path)
        ko_otu_table = 'ko_abundance.xls'
        enzyme_file_input = 'enzyme_abundance.xls'
        ko_outfile = os.path.join(out_dir, "prediction_KO.xls")
        module_outfile = os.path.join(out_dir, "prediction_module.xls")
        enzyme_outfile = os.path.join(out_dir, "prediction_enzyme.xls")
        level1_outfile = os.path.join(out_dir, "prediction_pathway.L1.xls")
        level2_outfile = os.path.join(out_dir, "prediction_pathway.L2.xls")
        level3_outfile = os.path.join(out_dir, "prediction_pathway.L3.xls")
        all_file_list = [ko_outfile, module_outfile, enzyme_outfile,level1_outfile, level2_outfile, level3_outfile]
        for file in all_file_list:
            if os.path.exists(file):
                os.remove(file)
        if otu_name.endswith(".gz"):
            otu_file = self.unzip_file(otu_path, ko_otu_table)
        else:
            otu_file = otu_table
        enzyme_path = enzyme
        enzyme_name = os.path.basename(enzyme_path)
        if enzyme_name.endswith(".gz"):
            enzyme_file_in = self.unzip_file(enzyme_path, enzyme_file_input)
        else:
            enzyme_file_in = enzyme
        # out_enzyme_file = "enzyme_file.xls"
        # new_enzyme_file_in = self.split_ec(enzyme_file_in, out_enzyme_file)
        data1 = pd.read_table(enzyme_file_in, sep="\t", header=0)
        sample_list1 = list(data1.columns)[1:]
        data1.columns = ["function"] + sample_list1
        data1['function'] = data1['function'].str.split(':', 1).str[1]
        ec_database_file = self.ec_database
        ec_database = pd.read_table(ec_database_file, sep='\t', header=0)
        ec_database["enzyme_id"] = ec_database["enzyme_id"].str.split(':', 1).str[1]
        merge_file1 = pd.merge(ec_database, data1, left_on='enzyme_id', right_on='function', how='inner')
        ###生成enzyme的丰度表  此表是picrust2库中注释得到的数据，本操作只是增加上了酶的描述信息
        enzyme_file_list = ["enzyme_id", "enzyme_category"] + sample_list1
        enzyme_file = merge_file1.loc[:, enzyme_file_list]
        enzyme_file.drop_duplicates("enzyme_id", "first", inplace=True)
        enzyme_file = enzyme_file[enzyme_file["enzyme_id"].astype('str') != "-"]
        enzyme_file.columns = ["Enzyme_ID", "Enzyme_Category"] + sample_list1
        enzyme_file.to_csv(enzyme_outfile, sep='\t', index=False, header=1)

        ##下面是通过找KO的对应关系找到对应的pathway水平，计算module的丰度
        data = pd.read_table(otu_file, sep="\t", header=0)
        sample_list = list(data.columns)[1:]
        data.columns = ["ko_id"] + sample_list
        database_file = self.database[database]
        data_database = pd.read_table(database_file, sep="\t", header=0)
        merge_file = pd.merge(data_database, data, left_on='ko_id', right_on='ko_id', how='inner')
        ###生成KO的丰度表
        ko_file_list = ["ko_id", "ko_desc"] + sample_list
        ko_file = merge_file.loc[:, ko_file_list]
        ko_file.drop_duplicates("ko_id", "first", inplace=True)
        ko_file = ko_file[ko_file["ko_id"] != "-"]
        ko_file.columns = ["KO_ID", "KO_Description"] + sample_list
        ko_file.to_csv(ko_outfile, sep='\t', index=False, header=1)
        ###生成module丰度表
        module_file_list = ["module_id"] + sample_list
        module_file = merge_file.loc[:, module_file_list]
        module = module_file['module_id'].str.split(';',expand=True).stack().to_frame()
        module = module.reset_index(level=1, drop=True).rename(columns={0:'module_id'})
        module_file = module_file.drop(['module_id'], axis=1).join(module)
        module_file = module_file.groupby('module_id').sum()
        module_file = module_file.reset_index()
        module_database = pd.read_table(self.module_database, sep="\t", header=0)
        new_module_file = pd.merge(module_file, module_database, left_on='module_id', right_on='ENTRY', how='inner')
        print(new_module_file.columns)
        module_file2_list = ["module_id", "NAME"] + sample_list
        new_module_file = new_module_file.loc[:, module_file2_list]
        new_module_file.drop_duplicates("module_id", "first", inplace=True)
        new_module_file = new_module_file[new_module_file["module_id"].astype('str') != "-"]
        new_module_file.columns = ["Module_ID", "Module_Category"] + sample_list
        new_module_file.to_csv(module_outfile, sep='\t', index=False, header=1)

        ##下面是通过生成小ko的id来找对应的level1、level2、level3层级对应关系
        pathway_file_list = ["pathway_id"] + sample_list
        pathway_file = merge_file.loc[:, pathway_file_list]
        pathway_all_file = pathway_file.drop("pathway_id", axis=1).join(pathway_file['pathway_id'].str.split(";", expand=True).stack().reset_index(level=1,drop=True).rename('pathway_id'))
        ko_level_file = self.ko_database
        ko_level_file_database = pd.read_table(ko_level_file, sep="\t", header=0)
        merge_level_file = pd.merge(ko_level_file_database, pathway_all_file, left_on='pathway', right_on='pathway_id', how='inner')

        ###生成pathway_L1的丰度表
        ko_file_list = ["level1"] + sample_list
        ko_level1_file = merge_level_file.loc[:, ko_file_list]
        # ko_level1_file.drop_duplicates("level1", "first", inplace=True)
        ko_level1_file = ko_level1_file.groupby("level1").sum()
        ko_level1_file = ko_level1_file.reset_index()
        ko_level1_file = ko_level1_file[ko_level1_file["level1"].astype('str') != "-"]
        ko_level1_file.columns = ["Level1"] + sample_list
        ko_level1_file.to_csv(level1_outfile, sep='\t', index=False, header=1)
        ###生成pathway_L2的丰度表
        ko_file_list = ["level2"] + sample_list
        ko_level2_file = merge_level_file.loc[:, ko_file_list]
        # ko_level2_file.drop_duplicates("level2", "first", inplace=True)
        ko_level2_file = ko_level2_file.groupby("level2").sum()
        ko_level2_file = ko_level2_file.reset_index()
        ko_level2_file = ko_level2_file[ko_level2_file["level2"].astype('str') != "-"]
        ko_level2_level1_file = pd.merge(ko_level2_file, ko_level_file_database, left_on='level2', right_on='level2', how='inner')
        print(ko_level2_level1_file.columns)
        ko_level2_level1_file_list = ["level2", "level1"] + sample_list
        ko_level2_level1_file = ko_level2_level1_file.loc[:, ko_level2_level1_file_list]
        ko_level2_level1_file.drop_duplicates("level2", "first", inplace=True)
        ko_level2_level1_file.columns = ["Level2", "Level1"] + sample_list
        ko_level2_level1_file.to_csv(level2_outfile, sep='\t', index=False, header=1)
        ###生成pathway_L3的丰度表
        ko_file_list = ["level3"] + sample_list
        ko_level3_file = merge_level_file.loc[:, ko_file_list]
        # ko_level3_file.drop_duplicates("level3", "first", inplace=True)
        ko_level3_file = ko_level3_file.groupby("level3").sum()
        ko_level3_file = ko_level3_file.reset_index()
        new_ko_level3_file = pd.merge(ko_level3_file, ko_level_file_database, left_on='level3', right_on='level3', how='inner')
        print(new_ko_level3_file.columns)
        ko2_file_list = ["level3","pathway", "level2", "level1"] + sample_list
        ko_level3_level2_level1_file = new_ko_level3_file.loc[:, ko2_file_list]
        ko_level3_file2 = ko_level3_level2_level1_file[ko_level3_level2_level1_file["level3"].astype('str') != "-"]
        ko_level3_file2.columns = ["level3", "Pathway_ID", "Level2", "Level1"] + sample_list
        ko_level3_file2.to_csv(level3_outfile, sep='\t', index=False, header=1)

    def run_metacyc_anno(self, otu_table, out_dir, database):
        """
        根据输入的otu表得到cog的注释结果
        :param otu_table: 输入的otu丰度表
        :param out_dir: 输出结果文件夹
        :param database: 输入数据库类型
        :return:
        """
        print("开始进行metacyc注释分析")
        otu_path = otu_table
        otu_name = os.path.basename(otu_path)
        metacyc_outfile = os.path.join(out_dir, "MetaCyc_pathway_pred.xls")
        metacyc_table = 'metacyc_abundance.xls'
        if os.path.exists(metacyc_outfile):
            os.remove(metacyc_outfile)
        if otu_name.endswith(".gz"):
            otu_file = self.unzip_file(otu_path, metacyc_table)
        else:
            otu_file = otu_table
        data = pd.read_table(otu_file, sep="\t", header=0)
        sample_list = list(data.columns)[1:]
        data.columns = ["pathway"] + sample_list
        database_file = self.database[database]
        data_database = pd.read_table(database_file, sep="\t", header=0)
        merge_file = pd.merge(data_database, data, left_on='pathway', right_on='pathway', how='inner')
        meta_file_list = ["pathway", "pathway_des"] + sample_list
        meta_file = merge_file.loc[:, meta_file_list]
        meta_file.drop_duplicates("pathway", "first", inplace=True)
        meta_file.columns = ["Pathway", "Pathway_Description"] + sample_list
        meta_file.to_csv(metacyc_outfile, sep='\t', index=False, header=1)

    def run_enzyme_anno(self, otu_table, out_dir, database):
        """
        对酶进行注释
        :param otu_table: 酶的丰度表
        :param out_dir: 输出文件夹
        :param database: 输入数据库类型
        :return:
        """
        print("开始进行kegg注释分析")
        otu_path = otu_table
        otu_name = os.path.basename(otu_path)
        enzyme_file_input = 'enzyme_abundance.xls'
        enzyme_outfile = os.path.join(out_dir, "prediction_enzyme.xls")
        if os.path.exists(enzyme_outfile):
            os.remove(enzyme_outfile)
        if otu_name.endswith(".gz"):
            enzyme_file_in = self.unzip_file(otu_path, enzyme_file_input)
        else:
            enzyme_file_in = otu_table
        data1 = pd.read_table(enzyme_file_in, sep="\t", header=0)
        sample_list1 = list(data1.columns)[1:]
        data1.columns = ["function"] + sample_list1
        data1['function'] = data1['function'].str.split(':', 1).str[1]
        ec_database_file = self.ec_database
        ec_database = pd.read_table(ec_database_file, sep='\t', header=0)
        ec_database["enzyme_id"] = ec_database["enzyme_id"].str.split(':', 1).str[1]
        merge_file1 = pd.merge(ec_database, data1, left_on='enzyme_id', right_on='function', how='inner')
        ###生成enzyme的丰度表  此表是picrust2库中注释得到的数据，本操作只是增加上了酶的描述信息
        enzyme_file_list = ["enzyme_id", "enzyme_category"] + sample_list1
        enzyme_file = merge_file1.loc[:, enzyme_file_list]
        enzyme_file.drop_duplicates("enzyme_id", "first", inplace=True)
        enzyme_file = enzyme_file[enzyme_file["enzyme_id"].astype('str') != "-"]
        enzyme_file.columns = ["Enzyme_ID", "Enzyme_Category"] + sample_list1
        enzyme_file.to_csv(enzyme_outfile, sep='\t', index=False, header=1)

    def unzip_file(self, file, out_file):
        """
        解压文件
        :param file: 待解压文件
        :param out_file: 解压后文件
        :return:
        """
        unzip_cmd = "zcat " +  file + " > " + out_file
        try:
            subprocess.check_output(unzip_cmd, shell=True)
            print("unzip done")
        except subprocess.CalledProcessError:
            raise Exception("unzip error")
        return out_file

    def split_ec(self,infile, outfile):
        """
        处理酶的丰度表
        :param infile: 酶丰度表
        :param outfile: 修改后的文件
        :return:
        """
        with open(infile, 'r') as f, open(outfile, 'w') as w:
            lines = f.readlines()
            head = lines[0]
            w.write(head)
            for line in lines[1:]:
                line = line.strip().split("\t")
                enzyme = line[0]
                new_enzyme = enzyme.split(":")
                if len(new_enzyme) != 2:
                    print(enzyme)
                n_enzyme = new_enzyme[1]
                new_line = line
                new_line[0] = n_enzyme
                abundance = "\t".join(new_line)
                w.write("{}\n".format(abundance))
        return outfile


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',metavar='[abundance_table]',required=True,help='Input abundance table ')
    parser.add_argument('-database',metavar='[database]',required=False,help='Input database type')
    parser.add_argument('-e',metavar='[enzyme]',required=False,help='Input enzyme type')
    parser.add_argument('-o',metavar='[output_dir]',required=True,help='output dir name')
    args = parser.parse_args()
    otu_table = args.i
    out_dir = args.o
    database = args.database
    meta_anno = PickAnnotation()
    if database == 'cog':
        meta_anno.run_cog_anno(otu_table, out_dir, database)
    elif database == 'metacyc':
        meta_anno.run_metacyc_anno(otu_table, out_dir, database)
    elif database == 'kegg':
        enzyme = args.e
        meta_anno.run_kegg_anno(otu_table, out_dir, database, enzyme)
    elif database == 'enzyme':
        enzyme = args.e
        meta_anno.run_enzyme_anno(otu_table, out_dir, database)