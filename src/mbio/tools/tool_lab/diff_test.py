# -*- coding: utf-8 -*-


import os, re, subprocess
from biocluster.agent import Agent
from biocluster.tool import Tool
import pandas as pd
import numpy as np
from biocluster.core.exceptions import OptionError
from mbio.packages.metabolome.scripts.metastat import two_group_test
from mbio.packages.metabolome.scripts.merge_table import merge_table
#from mbio.packages.tool_lab.get_box import get_box_value


class DiffTestAgent(Agent):
    """
    两组差异检验，student T检验，welch T检验，wilcox秩和检验
    """

    def __init__(self, parent):
        super(DiffTestAgent, self).__init__(parent)
        options = [
            {'name': 'exp_file', 'type': 'infile', 'format': 'metabolome.express,metabolome.metab_abun'},  # 表达矩阵文件
            {'name': 'test_method', 'type': 'string', 'default': 't-test'},  # 差异检验方法
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},  # 分组文件
            {'name': 'group_name', 'type': 'string', 'default': ''},  # 用于检测的组名，eg："A|B" ,两两比较
            {'name': 'side_type', 'type': 'string', 'default': 'two-tailed'},
            # 单尾或双尾检验 two-tailed,left-tailed,right-tailed
            {'name': 'ci', 'type': 'float', 'default': 0.95},  # 显著性检验水平
            {'name': 'mul_test', 'type': 'string', 'default': 'none'}
            # ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"]
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("exp_file").is_set:
            raise OptionError("请传入丰度矩阵！")
        if not self.option("group_file").is_set:
            raise OptionError("请传入group文件！")
        if not self.option("test_method") in ["t-test", 'welch', 'wilcox']:
            raise OptionError("仅支持t-test, welch, wilcox两组检验方法！")
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '5G'

    def end(self):
        super(DiffTestAgent, self).end()


class DiffTestTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(DiffTestTool, self).__init__(config)
        self.r_path = "/program/R-3.3.1/bin/Rscript"   #change to R-3.3.1




    def run(self):
        """
        运行
        """
        super(DiffTestTool, self).run()
        self.logger.info("开始运行命令！")
        self.logger.info(self.option("exp_file").prop["path"])

        self.diff_g_path, self.target_group_info = self.make_group(self.option("group_file").prop["path"], self.option("group_name"))

        self.run_diff()
        self.treat_result()
        self.set_output()

    def make_group(self, group_file, diff_group):
        """
        根据差异分组名，从总group表中生成各差异组group
        :param group_file: 总group文件
        """

        group_dict = {}
        with open(group_file, "r") as f:
            for line in f:
                line = line.strip()
                line1 = line.split("\t")
                if not "#" in line:
                    group = line1[1]
                    name = line1[0]
                    if not group_dict.has_key(group):
                        group_dict[group] = []
                    group_dict[group].append(name)


        two_groups = diff_group.split("|")
        two_groups.reverse()

        diff_g_path = self.work_dir + "/new_group.xls"

        target_group_info = {}
        with open(diff_g_path, "w") as gf:
            gf.write("#sample\tgroup\n")
            for i in two_groups:
                if not group_dict.has_key(i):
                    raise OptionError('分组样本不在丰度表所拥有的样本内，请检查分组方案')
                names = group_dict[i]
                target_group_info[i] = names
                for j in names:
                    gf.write(j + "\t" + i + "\n")
        return diff_g_path,target_group_info

    def run_diff(self):
        """
        两组差异分析
        """

        if self.option("test_method") == "t-test":
            test_method = "student"
        elif self.option("test_method") == "welch":
            test_method = "welch"
        elif self.option("test_method") == "wilcox":
            test_method = "wilcox"
        else:
            test_method = "student"
        if self.option("side_type") == "two-tailed":
            side_type = "two.sided"
        elif self.option("side_type") == "left-tailed":
            side_type = "less"
        elif self.option("side_type") == "right-tailed":
            side_type = "greater"
        else:
            side_type = "two.sided"

        myfile = self.diff_g_path
        group = myfile.split("/")[-1].split("_group")[0]
        exp = self.option("exp_file").prop["path"]
        self.check_group(myfile)
        out = self.work_dir + "/result.xls"
        box_out =  self.work_dir + "/boxfile.xls"
        ci = self.option("ci")
        two_group_test(exp, myfile, out, box_out, test_method, ci=1 - ci,
                       test_type=side_type, mul_test=self.option("mul_test"), norm="F")
        os.rename(self.work_dir + "/run_{}_test.r".format(test_method), self.work_dir +
                  "/run_{}_test".format(test_method) + group + ".r")
        cmd = self.r_path + " run_{}_test".format(test_method) + group + ".r"


        self.result = out
        self.group_file = myfile


        command = self.add_command("diff" , cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run {} succeed".format(cmd))
        else:
            self.set_error("run %s failed", variables=(cmd))
            raise Exception("run {} failed".format(cmd))





    def treat_result(self):
        exp_data = pd.read_table(self.option('exp_file').path, sep='\t',index_col=0)
        final_data = []
        for g in  self.target_group_info:
            tmp_data = exp_data[self.target_group_info[g]]
            box = tmp_data.apply(lambda x: x.describe(),axis=1)
            delta_q = box['75%'] - box['25%']
            tmp_data['up'] = box['75%'] + delta_q*1.5
            tmp_data['down'] = box['25%'] - delta_q*1.5

            sub_min = []
            sub_max = []
            abnormal = []
            ##找异常点
            for index in tmp_data.index:
                up = tmp_data['up'][index]
                down = tmp_data['down'][index]
                row = tmp_data.loc[index]
                ab = row[row.map(lambda x: True if x>up or x<down else False)]
                nor = row[row.map(lambda x: False if x>up or x<down else True)]
                sub_min.append(nor.min())
                sub_max.append(nor.max())
                abnormal.append(';'.join([str(s[0])+':'+str(s[1]) for s in zip(ab.index.tolist(),ab.tolist())]))

            box.drop(['count'],axis=1,inplace=True)
            box[g+'--sub-min'] = sub_min
            box[g+'--sub-max'] = sub_max
            box[g+'--abnormal'] = abnormal
            box.rename(columns={'mean':g+'--mean','std':g+'--std','min':g+'--min','max':g+'--max','25%':g+'--25%','75':g+'--75%', '50%':g+'--50%'},inplace=True)
            final_data.append(box)

        final = pd.concat(final_data,axis=1)

        two_groups = self.option('group_name').split('|')
        def my_test(row):
            if row[0] != 0:
                fc = float(row[1]) / float(row[0])
            else:
                if row[1] != 0:
                    fc = ''
                else:
                    fc = 0
            return fc

        final["FC"] = final[[two_groups[0]+'--mean', two_groups[1]+'--mean']].apply(my_test, axis=1)
        final.to_csv('box_result.xls',sep='\t')


    def keep_group_samples(self, group_file):
        group_samples = []
        with open(group_file, "r") as f:
            for line in f:
                line = line.strip().split("\t")
                if not "#" in line[0]:
                    sample = line[0]
                    if not sample in group_samples:
                        group_samples.append(sample)
        return group_samples

    def set_output(self):
        """

        """
        self.logger.info("set output")
        self.end()

    def check_group(self, groupfile):
        with open(groupfile, 'rb') as g:
            info = g.readlines()
            sample = []
            group_dict = {}
            group_num = {}
            for line in info[1:]:
                sample.append(line.strip('\n').split('\t')[0])
                line = line.split()
                group_dict[line[0]] = line[1]
            gnames = group_dict.values()
            for gname in gnames:
                group_num[gname] = gnames.count(gname)
            sam_num = group_num.values()
            if 1 in sam_num:
                self.set_error('一个分组类别下面至少有两个样本，请重新分组')
            if len(sample) != len(set(sample)):
                self.set_error('不同的分组类别下有相同的样本，此分析不支持该分组方案，请重新分组')


