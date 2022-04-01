# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2018.06.05

import os, re, subprocess
from biocluster.agent import Agent
from biocluster.tool import Tool
import pandas as pd
import numpy as np
from biocluster.core.exceptions import OptionError
from mbio.packages.metabolome.scripts.metastat import two_group_test

import traceback
from mbio.packages.metabolome.scripts.merge_table import merge_table
from mbio.packages.metabolome.common import Relation


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
            {'name': 'group_name', 'type': 'string', 'default': ''},  # 用于检测的组名，eg："A|B";"A|C" ,两两比较
            {'name': 'side_type', 'type': 'string', 'default': 'two-tailed'},
            # 单尾或双尾检验 two-tailed,left-tailed,right-tailed
            {'name': 'ci', 'type': 'float', 'default': 0},  # 显著性检验水平
            {'name': 'mul_test', 'type': 'string', 'default': 'none'}
            # ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"]
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("exp_file").is_set:
            raise OptionError("请传入丰度矩阵！", code="34700301")
        if not self.option("group_file").is_set:
            raise OptionError("请传入group文件！", code="34700302")
        if not self.option("test_method") in ["t-test", 'welch', 'wilcox', "t-test-paired", 'welch-paired', 'wilcox-paired']:
            raise OptionError("仅支持t-test, welch, wilcox两组检验方法！", code="34700303")
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
        #self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.3/bin/Rscript"
        self.r_path = "/program/R-3.3.1/bin/Rscript"   #change to R-3.3.1
        self.cmds = []
        self.results = {}
        self.group_file = {}

    def run(self):
        """
        运行
        """
        super(DiffTestTool, self).run()
        self.logger.info("开始运行命令！")
        self.logger.info(self.option("exp_file").prop["path"])
        if self.option("group_name") != "":
            self.diff_g_paths = self.make_group(self.option("group_file").prop["path"], self.option("group_name"))
        else:
            self.diff_g_paths = [self.option("group_file").prop["path"]]
        self.logger.info(self.diff_g_paths)
        self.run_diff()
        self.treat_result()
        self.set_output()

    def make_group(self, group_file, diff_group):
        """
        根据差异分组名，从总group表中生成各差异组group
        :param group_file: 总group文件
        :param diff_group: 差异组名， eg："A|B";"A|C"
        :return: 返回生成差异group文件路径diff_g_paths列表
        """
        diff_groups = diff_group.split(";")
        group_dict = {}
        diff_g_paths = []
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
        for each in diff_groups:
            #two_groups = each.split("|")
            #diff_name = each.replace("|", "_vs_")
            two_groups = each.split("_vs_")
            two_groups.reverse()    # modify 20181106 other/control
            diff_name = each
            diff_g_path = self.work_dir + "/" + diff_name + "_group"
            diff_g_paths.append(diff_g_path)
            with open(diff_g_path, "w") as gf:
                gf.write("#sample\tgroup\n")
                for i in two_groups:
                    if not group_dict.has_key(i):
                        raise OptionError('分组样本不在丰度表所拥有的样本内，请检查分组方案', code="34700304")
                    names = group_dict[i]
                    for j in names:
                        gf.write(j + "\t" + i + "\n")
            self.logger.info("生成{}_group表成功".format(each))
        return diff_g_paths

    def run_diff(self):
        """
        两组差异分析
        """
        count = 0
        if self.option("test_method") == "t-test":
            test_method = "student"
        elif self.option("test_method") == "welch":
            test_method = "welch"
        elif self.option("test_method") == "wilcox":
            test_method = "wilcox"
        elif self.option("test_method") == "t-test-paired":
            test_method = "student-paired"
        elif self.option("test_method") == "welch-paired":
            test_method = "welch-paired"
        elif self.option("test_method") == "wilcox-paired":
            test_method = "wilcox-paired"

        if self.option("side_type") == "two-tailed":
            side_type = "two.sided"
        elif self.option("side_type") == "left-tailed":
            side_type = "less"
        elif self.option("side_type") == "right-tailed":
            side_type = "greater"
        for myfile in self.diff_g_paths:
            group = myfile.split("/")[-1].rpartition("_group")[0]
            exp = self.option("exp_file").prop["path"]
            #exp = self.remove_zero(exp, myfile)
            self.check_group(myfile)
            out = self.work_dir + "/" + group + '_result' + '.xls'
            box_out =  self.work_dir + "/" + group + '_boxfile.xls'
            ci = self.option("ci")
            two_group_test(exp, myfile, out, box_out, test_method, ci=1 - ci,
                           test_type=side_type, mul_test=self.option("mul_test"), norm="F")
            os.rename(self.work_dir + "/run_{}_test.r".format(test_method), self.work_dir +
                      "/run_{}_test".format(test_method) + group + ".r")
            cmd = self.r_path + " run_{}_test".format(test_method) + group + ".r"
            self.cmds.append(cmd)
            self.results[group] = out
            self.group_file[group] = myfile
            count += 1
            command = self.add_command("diff" + str(count), cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("run {} succeed".format(cmd))
            else:
                self.set_error("run %s failed", variables=(cmd), code="34700301")
                raise Exception("run {} failed".format(cmd))
        if count == len(self.cmds):
            self.logger.info("全部命令运行完成！！！")
        else:
            self.set_error("cmd未全部完成", code="34700303")
            raise Exception("cmd未全部完成")

    def treat_result(self):
        """
        处理文件为所需输出格式，head: metab\tFC\tP_value
        """
        exp_file = self.option("exp_file").prop["path"]
        exp_table = pd.read_table(exp_file, sep="\t", header=0)
        def my_test(row):
            if row[1] != 0:
                fc = float(row[2]) / float(row[1])
            else:
                if row[2] != 0:
                    fc = None
                else:
                    fc = 0
            return fc

        for group in self.results.keys():
            each = self.results[group]
            group_sample = self.keep_group_samples(self.group_file[group])
            group_sample = ["metab_id"] + group_sample
            group_exp = exp_table[group_sample]
            group_exp.rename(columns={group_exp.columns[0]: "Metab"}, inplace=True)
            if not os.path.exists(each):
                self.set_error("%s未正确生成!", variables=(each), code="34700305")
                raise Exception("{}未正确生成!".format(each))
            df_tmp = pd.read_table(each, sep="\t", header=0)
            # df.drop(df.columns[[2, 4, 6]], axis=1, inplace=True) # test 结果表列变化相应修改
            mean_and_sd = df_tmp[df_tmp.columns[0:5]]
            df_tmp.drop(df_tmp.columns[[2, 4]], axis=1, inplace=True)
            df = df_tmp.loc[:,df_tmp.columns[[0,1,2]]]
            df["pvalue"] = df_tmp["pvalue"]
            df['fdr'] = df_tmp["corrected_pvalue"]  ##
            df.rename(columns={df.columns[0]: "Metab", df.columns[3]: "P_value"}, inplace=True)
            df["FC"] = df.apply(my_test, axis=1)
            df = df.drop(df.columns[[1, 2]], axis=1)
            fc_max = 2 * (float(df["FC"].max()))  ## FC分子为0时设为最大FC的2倍
            df["FC"] = df["FC"].replace(np.nan, fc_max)
            outfile = each.replace(self.work_dir, self.output_dir)
            merge_table(df, group_exp, "Metab", outfile=outfile, how="inner")
            #df.to_csv(outfile, sep='\t', index=False)

            ## produce  plot box data and bar std data
            box = pd.read_table(self.work_dir+'/'+group+'_boxfile.xls',sep='\t',header=0)
            box.rename(columns={box.columns[0]:"Metab"},inplace=True)
            mean_and_sd.rename(columns={mean_and_sd.columns[0]:"Metab"},inplace=True)
            plot_data = pd.merge(mean_and_sd,box,on='Metab')
            if not os.path.exists(self.work_dir+'/plot_dir'):
                os.mkdir(self.work_dir+'/plot_dir')
            plot_data.to_csv(self.work_dir+'/plot_dir/'+group+'_plot.xls',sep='\t',index=False)

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
        id 转化
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
            # if 1 in sam_num:
            #     self.set_error('一个分组类别下面至少有两个样本，请重新分组', code="34700305")
            if len(sample) != len(set(sample)):
                self.set_error('不同的分组类别下有相同的样本，此分析不支持该分组方案，请重新分组', code="34700306")

    def trans_data(self, oldfile, newfile, map_table, id_name_dict):
        self.metab_trans.get_trans_file(map_table, id_name_dict, oldfile, newfile)
