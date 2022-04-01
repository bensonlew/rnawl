#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# @Author  : houshuang 2019/9/27

import os
# import subprocess
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import shutil
import subprocess
import re


class Tax4funAgent(Agent):
    def __init__(self, parent):
        super(Tax4funAgent, self).__init__(parent)
        options = [
            {"name": "database", "type": "string"},
            {"name": "in_otu_table", "type": "infile", "format": "meta.otu.otu_table"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数检测
        """
        if not os.path.isfile(self.option("in_otu_table").prop['path']):
            raise OptionError("ASV表不存在!")

    def set_resource(self):
        """
        设置所需资源
        """
        self.cpu = 1
        self.memory = '5G'

    def end(self):
        super(Tax4funAgent, self).end()


class Tax4funTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(Tax4funTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.database = self.config.SOFTWARE_DIR + '/database/meta/tax4fun/SILVA123'
        database_dict = {'silva123/16s_bacteria': ['SILVA123_16s_bacteria_ref_uniq.tax', 'SILVA123_16s_bacteria_origi_uniq.tax'],
                         'silva123/16s_archaea': ['SILVA123_16s_archaea_ref_uniq.tax', 'SILVA123_16s_archaea_origi_uniq.tax'],
                         'silva123/16s': ['SILVA123_16s_ref_uniq.tax', 'SILVA123_16s_origi_uniq.tax'],
                         'silva119/16s_bacteria': ['SILVA119_16s_bacteria_ref_uniq.tax', 'SILVA119_16s_bacteria_origi_uniq.tax'],
                         'silva119/16s_archaea': ['SILVA119_16s_archaea_ref_uniq.tax', 'SILVA119_16s_archaea_origi_uniq.tax'],
                         'silva119/16s': ['SILVA119_16s_ref_uniq.tax', 'SILVA119_16s_origi_uniq.tax'],
                         'silva128/16s_archaea': ['SILVA128_16s_archaea_ref_uniq.tax', 'SILVA128_16s_archaea_origi_uniq.tax'],
                         'silva128/16s_bacteria': ['SILVA128_16s_bacteria_ref_uniq.tax', 'SILVA128_16s_bacteria_origi_uniq.tax'],
                         'silva128/16s': ['SILVA128_16s_ref_uniq.tax', 'SILVA128_16s_origi_uniq.tax'],
                         'silva132/16s_archaea': ['SILVA132_16s_archaea_ref_uniq.tax', 'SILVA132_16s_archaea_origi_uniq.tax'],
                         'silva132/16s_bacteria': ['SILVA132_16s_bacteria_ref_uniq.tax', 'SILVA132_16s_bacteria_origi_uniq.tax'],
                         'silva132/16s': ['SILVA132_16s_ref_uniq.tax', 'SILVA132_16s_origi_uniq.tax'],
                         'silva138/16s_archaea': ['SILVA138_16s_archaea_ref_uniq.tax','SILVA138_16s_archaea_origi_uniq.tax'],
                         'silva138/16s_bacteria': ['SILVA138_16s_bacteria_ref_uniq.tax','SILVA138_16s_bacteria_origi_uniq.tax'],
                         'silva138/16s': ['SILVA138_16s_ref_uniq.tax', 'SILVA138_16s_origi_uniq.tax'],
                         'nt_v20200327/16s_archaea': ['nt_v20200327_16s_archaea_ref_uniq.tax', 'nt_v20200327_16s_archaea_origi_uniq.tax'],
                         'nt_v20200327/16s_bacteria': ['nt_v20200327_16s_bacteria_ref_uniq.tax', 'nt_v20200327_16s_bacteria_origi_uniq.tax'],
                         'nt_v20200327/16s': ['nt_v20200327_16s_ref_uniq.tax', 'nt_v20200327_16s_origi_uniq.tax'],
                         'greengenes135/16s': ['greengene_16s_ref_uniq.tax', 'greengene_16s_origi_uniq.tax'],
                         'greengenes135/16s_archaea': ['greengene_16s_archaea_ref_uniq.tax', 'greengene_16s_archaea_origi_uniq.tax'],
                         'greengenes135/16s_bacteria': ['greengene_16s_bacteria_ref_uniq.tax', 'greengene_16s_bacteria_origi_uniq.tax'],
                         'rdp11.5/16s': ['rdp_16s_ref_uniq.tax', 'rdp_16s_origi_uniq.tax'],
                         'rdp11.5/16s_bacteria': ['rdp_16s_bacteria_ref_uniq.tax', 'rdp_16s_bacteria_origi_uniq.tax'],
                         'rdp11.5/16s_archaea': ['rdp_16s_archaea_ref_uniq.tax', 'rdp_16s_archaea_origi_uniq.tax']}

        if self.option('database') in ['silva119/16s_bacteria', 'silva119/16s_archaea', 'silva119/16s']:
            self.database = self.config.SOFTWARE_DIR + '/database/meta/tax4fun/SILVA119'
        else:
            self.database = self.config.SOFTWARE_DIR + '/database/meta/tax4fun/SILVA123'
        self.modified_taxon = self.config.SOFTWARE_DIR + '/database/taxon_db/origi/uniq_name/' + database_dict[
            str(self.option('database'))][0]
        self.origi_taxon = self.config.SOFTWARE_DIR + '/database/taxon_db/origi/uniq_name/' + database_dict[
            str(self.option('database'))][1]

        self.R_path = '/program/R-3.3.1/bin/'
        self.R_path_2 = software_dir + '/program/R-3.3.1/bin/'
        self.tax4fun_r = self.config.PACKAGE_DIR + '/meta/scripts/tax4fun.R'
        self.outfile = self.work_dir + '/predictions_KO.xls'
        self.kegg = software_dir + '/database/meta/tax4fun/kegg_v94.2_pathway.anno.txt'
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
        # self.enzyme_database = software_dir + '/database/KEGG/kegg_2017-05-01/files/enzyme.txt'
        self.enzyme_database = software_dir + "/database/meta/tax4fun/kegg_v94.2_enzyme.txt" ##更新数据库 20201214

    def run_tax4fun(self):
        """
        调用脚本
        """
        infile = os.path.join(self.work_dir, "out.txt")
        cmd = self.R_path_2 + 'Rscript {} -i {} -o {} -d {} -k {}'.format(self.tax4fun_r, infile, self.outfile,
                                                                        self.database, self.kegg)
        self.logger.info("运行命令：".format(cmd))
        self.logger.info(cmd)
        # command = self.add_command("tax4fun_cmd", cmd)
        # command.run()
        # self.wait(command)
        # if command.return_code == 0:
        #       subprocess.check_output(cmd, shell=True)
        #     self.logger.info("运行tax4fun_cmd正常完成")

        try:
            self.logger.info(cmd)
            subprocess.check_output(cmd, shell=True)
            self.logger.info("运行tax4fun_cmd正常完成")
        except:
            outfile = os.path.join(self.work_dir, "new_out.txt")
            new_infile = self.replace_0_otu(infile, outfile)
            cmd2 = self.R_path_2 + 'Rscript {} -i {} -o {} -d {} -k {}'.format(self.tax4fun_r, new_infile, self.outfile, self.database, self.kegg)
            self.logger.info(cmd2)
            try:
                subprocess.check_output(cmd2, shell=True)
                self.logger.info("第二次运行tax4fun完成")
            except Exception,e:
                self.set_error("运行tax4fun_cmd出错!",code="32707201")

    def replace_0_otu(self, infile, outfile):
        """
        1.首先过滤掉out.txt丰度表中所有丰度加和为0的otu；
        2.将otu表中的0进行替换,替换为全表最小值的万分之一，如果没有就直接替换为万分之一;
        3.进行file检查是否是有多行；
        """
        second_list = []

        with open(infile, 'r') as f, open(outfile, 'w') as w:
            lines = f.readlines()
            first_line = lines[0]
            w.write(first_line)
            all_name = lines[1].strip().split("\t")
            sample_name = all_name[1:-1]
            second_line = lines[1]
            w.write(second_line)
            line_number = 0
            for line in lines[2:]:##求最小值
                line = line.strip().split('\t')
                min_list = []
                for i in range(1, len(line)-1):
                    if float(line[i]) != 0.0:
                        line_min = float(line[i])
                        min_list.append(line_min)
                if len(min_list) != 0:
                    min_line = min(min_list)
                    second_list.append(min_line)
                else:
                    pass
                if not re.search(r"No blast hit", line[-1]):
                    line_number += 1
            self.logger.info(line_number)
            if line_number < 2:
                self.set_error("所注释到的物种结果少于2条不能进行Tax4fun分析")

            if len(second_list) != 0:
                min_table = float(min(second_list))##得到全表最小值
            else:
                min_table = float(1.0)
            for line in lines[2:]:
                line = line.strip().split('\t')
                data_list = []
                sum_all_data = 0
                for i in range(1, len(line)-1):
                    if float(line[i]) == 0.0:
                        line[i] = float(line[i]) + float(min_table) / 1000000
                    else:
                        line[i] = float(line[i])
                    sum_all_data += float(line[i])
                    data_list.append(line[i])
                if sum_all_data != 0.0:
                    w.write("{}\t{}\t{}\n".format("\t".join(line[0:1]), "\t".join(str(i) for i in data_list), line[-1]))
                else:
                    pass
        return outfile

    def generate_file(self):
        input_file = os.path.join(self.work_dir, "predictions_KO.xls.pathway.anno")
        pathway_file = pd.read_csv(input_file, sep='\t', encoding="utf-8")
        cols = list(pathway_file)
        cols.insert(2, cols.pop(cols.index('Pathway_level3')))
        data = pathway_file.loc[:, cols]
        data.to_csv("predictions_ko.L3.xls", header=True, index=False, sep="\t")

        data.drop("Level3_description", axis=1, inplace=True)
        data.drop("Pathway_level3", axis=1, inplace=True)
        group = data.groupby("Pathway_level2").sum()
        group = group.reset_index()
        new_data = data.loc[:, ["Pathway_level1", "Pathway_level2"]]
        new_data.drop_duplicates(inplace=True)
        out = pd.merge(new_data, group)
        out.to_csv("predictions_ko.L2.xls", header=True, index=False, sep="\t")
        # shutil.copy2("predictions_ko.L2.txt", "predictions_ko.L2.xls")

        data.drop("Pathway_level2", axis=1, inplace=True)
        group2 = data.groupby("Pathway_level1").sum()
        group2 = group2.reset_index()
        group2.to_csv("predictions_ko.L1.xls", header=True, index=False, sep="\t")

        ko = open(self.outfile, 'r')
        enzyme = open("predictions_Enzyme.xls", 'w')
        lines = ko.readlines()
        title = re.split("\t", lines[0])
        enzyme.write("Enzyme\tDescription\t" + "\t".join(title[2:]))
        ## 下面代码是为做酶的丰度的迭代和累加
        enzyme_description = {}
        with open(self.enzyme_database, "r") as enzy:
            enzy.readline()
            for line in enzy:
                line = line.strip().split("\t")
                enzyme_id = line[0]
                if enzyme_id not in enzyme_description:
                    enzyme_description[enzyme_id] = line[1:]

        enzyme_list = []
        enzyme_dict = {}
        for line in lines[1:]:
            if re.search("\[EC:", line):
                split = re.split("\t", line.strip())
                temp = re.split("\[EC:", split[1])
                if temp[0] == " ":
                    temp[0] = "-"  # decription不存在时用-代替
                temp[1] = temp[1].rstrip('\"') ## 酶正则匹配含有特殊符号双引号
                temp[1] = temp[1].rstrip('\]') ## 酶正则匹配含有特殊符号中括号
                temp[1] = re.sub(" ", ";", temp[1])  # 酶的编号有两个时，用";"分隔
                # all_temp = temp[1] + "\t" + temp[0].strip()
                temp_list = temp[1].split(";")
                if len(temp_list) > 0:
                    for enz in temp_list:
                        if enz not in enzyme_list:
                            enzyme_list.append(enz)
                            enzyme_dict[enz] = "\t".join(split[2:])
                        else:
                            if enz in enzyme_dict:
                                origin_abundance = enzyme_dict[enz].strip().split("\t")
                                now_list = split[2:]
                                new_list = [float(origin_abundance[i])+float(now_list[i]) for i in range(len(origin_abundance))]
                                new_list = [str(j) for j in new_list]
                                enzyme_dict[enz] = "\t".join(new_list)


        for key in enzyme_list:
            if key in enzyme_dict:
                if key in enzyme_description:
                    des = enzyme_description[key][-1].strip()
                    enzyme.write("{}\t{}\t{}\n".format(key, des,enzyme_dict[key].strip()))
                # else: ## 产品线需求如果酶的层级不为4个不参与统计
                #     key_list = key.split("\.")
                #     n = 0
                #     index_number = 0
                #     for key_site in key_list:
                #         if key_site == "-":
                #             n += 1
                #             if n == 1:
                #                 index_number = key_list.index(key_site)
                #     try:
                #         new_key = key.replace("-", 1)
                #     except:
                #         new_key = key
                #     if new_key in enzyme_description:
                #         if index_number != 0:
                #             des = enzyme_description[enzyme_description][index_number-1].strip()
                #     else:
                #         des = "-"
                #     enzyme.write("{}\t{}\t{}\n".format(key, des,enzyme_dict[key].strip()))
        ko.close()
        enzyme.close()

    def restore_spe_name(self):
        modify = open(self.modified_taxon, 'r')
        lines = modify.readlines()
        modify.close()
        modify_info = {}
        for line in lines:
            line = re.split("\t", line.strip())
            modify_info[line[1]] = line[0]
        origi = open(self.origi_taxon, 'r')
        lines2 = origi.readlines()
        origi.close()
        origi_info = {}
        for line in lines2:
            line = re.split("\t", line.strip())
            origi_info[line[0]] = line[1]
        f = open(self.option('in_otu_table').prop['path'], 'r')
        lines = f.readlines()
        f.close()
        out = open("out.txt", 'w')
        record = open("record.txt", 'w')
        out.write("# Constructed from biom file\n")
        out.write("# " + lines[0].strip() + "\ttaxonomy\n")
        count = 0
        for line in lines[1:]:
            count += 1
            split = re.split("\t", line.strip())
            temp = split[0]
            temp = re.sub("; ", ";", temp)
            temp = re.sub(";ASV.*", "", temp)
            if temp not in modify_info.keys():
                record.write(str(count) + ":" + temp + "\n")
                temp = "No blast hit"
            elif modify_info[temp] not in origi_info.keys():
                record.write("not found in origi: " + temp + "\n")
                temp = "No blast hit"
            else:
                temp = origi_info[modify_info[temp]]
            split = "\t".join(split[1:])
            out.write("ASV" + str(count) + "\t" + split + "\t" + temp + ";\n")
        out.close()
        record.close()

    def set_output(self):
        """
        将结果文件链接至output
        """
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        os.link(self.outfile, self.output_dir + '/predictions_KO.xls')
        os.link(self.work_dir + '/predictions_Enzyme.xls', self.output_dir + '/predictions_Enzyme.xls')
        os.link(self.work_dir + '/predictions_ko.L1.xls', self.output_dir + '/predictions_ko.L1.xls')
        os.link(self.work_dir + '/predictions_ko.L2.xls', self.output_dir + '/predictions_ko.L2.xls')
        os.link(self.work_dir + '/predictions_ko.L3.xls', self.output_dir + '/predictions_ko.L3.xls')
        self.logger.info("done")

    def run(self):
        """
        运行
        """
        super(Tax4funTool, self).run()
        self.logger.info("开始准备文件！")
        self.restore_spe_name()
        self.logger.info("开始运行tax4fun！")
        self.run_tax4fun()
        self.logger.info("生成其他文件！")
        self.generate_file()
        self.logger.info("set output")
        self.set_output()
        self.end()