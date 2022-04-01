# !usr/bin/python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.tool_lab.metastat import *
from mbio.packages.statistical.twogroup_CI import *
from mbio.packages.statistical.twosample_CI import *
from mbio.packages.statistical.mul_posthoc import *
from mbio.packages.tool_lab.get_box import get_box_value_row
import pandas as pd
import os
import re
import glob
from mbio.packages.metagenomic.common import check_command


class MannTestAgent(Agent):
    """
    statistical metastat+ 调用metastat.py 包进行差异性分析
    """
    def __init__(self, parent):
        super(MannTestAgent, self).__init__(parent)
        options = [
            {"name": "mann_input", "type": "infile", "format": "tool_lab.table"},  # 秩和检验的输入文件
            {"name": "mann_group", "type": "infile", "format": "tool_lab.group_table"},  # 秩和检验的输入分组文件
            {"name": "sample", "type": "string", "default": "column"},  # 样本名为列标签或行标签
            {"name": "compare_type", "type": "string", "default": "multi"},  # 比较策略
            {"name": "group_name1", "type": "string"},  # 秩和检验样本组1名称
            {"name": "group_name2", "type": "string"},  # 秩和检验样本组2名称
            {"name": "group_name", "type": "string"},  # 秩和检验样本组名称
            {"name": "mann_ci", "type": "float", "default": 0.05},  # 秩和检验的显著性水平
            {"name": "correction", "type": "string", "default": "none"},  # 秩和检验的多重检验校正
            {"name": "mann_type", "type": "string", "default": "two.side"},  # 秩和检验的选择单尾或双尾检验
            {"name": "kru_methor", "type": "string", "default": "tukeykramer"},  # 多组检验方法
            {"name": "coverage", "type": "float", "default": 0.95},  # 秩和检验的置信度
            {"name": "method", "type": "string", "default": "T"},
        ]
        self.add_option(options)

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
                raise OptionError('一个分组类别下面至少有两个样本，请重新分组')
            if len(sample) != len(set(sample)):
                raise OptionError('不同的分组类别下有相同的样本，此分析不支持该分组方案，请重新分组')
            gnames = list(set(gnames))
            if len(gnames) < 2:
                raise OptionError('至少有两组样本')

    def check_options(self):
        """
        检查参数设置
        """
        if not self.option("mann_input").is_set:
            raise OptionError('必须设置秩和检验输入的文件')
        if not self.option("mann_group").is_set:
            raise OptionError('必须设置秩和检验输入的分组文件')
        if self.option("correction") not in ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"]:
            raise OptionError('该多重检验校正的方法不被支持')
        if self.option("kru_methor") not in ["tukeykramer", "gameshowell", "welchuncorrected", "scheffe","dunn", "nemenyi", "conover-iman", "steel_dwass"]:
            raise OptionError('该多重检验校正的方法不被支持')
        if self.option("mann_ci") <= 0 or self.option("mann_ci") >= 1:
            raise OptionError('所输入的显著水平不在范围值内')
        if self.option("mann_type") not in ["two.side", "greater", "less"]:
            raise OptionError('所输入的类型不在范围值内')
        if self.option("coverage") not in [0.90, 0.95, 0.98, 0.99, 0.999]:
            raise OptionError('秩和检验的置信区间的置信度不在范围值内')
        if self.option("compare_type") == "multi":
            if self.option("group_name").find(";") >= 0:
                g_list = self.option("group_name").split(';')
            elif self.option("group_name").find(",") >= 0:
                g_list = self.option("group_name").split(',')
            else:
                raise OptionError('分组名请按英文字符“;”或“,”分隔！')
            g_list = list(set(g_list))
            new_list = [i for i in g_list if i != '']
            if len(new_list) < 3:
                raise OptionError('至少选择三个不同分组进行比较')

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            [r".*_result\.xls", "xls", "物种组间差异显著性比较结果表"],
            [r".*_boxfile\.xls", "xls", "组间差异显著性比较用于画箱线图的数据，包含四分位值"],
            [r".*_CI\.xls", "xls", "组间差异显著性比较两组，两样本比较的置信区间值以及效果量"],
            [r".*_bar\.xls", "xls", "单物种柱状图数据"],
            [r".*_kru_H\.xls", "xls", "多组间两两比较数据"]
        ])
        super(MannTestAgent, self).end()


class MannTestTool(Tool):
    def __init__(self, config):
        super(MannTestTool, self).__init__(config)
        self._version = "v1.0.1"
        self.package_path = 'packages/metastat.py'
        self.r_path = '/program/R-3.3.1/bin/Rscript'

    def run(self):
        """
        运行
        """
        super(MannTestTool, self).run()
        self.trans_sample()
        self.get_sub_group()
        if self.option('compare_type') != "multi":
            self.run_mann()
        else:
            self.run_kru()
        self.get_box()
        self.set_output()
        self.end()

    def trans_sample(self):
        """
        转置
        """
        sample = self.option("mann_input").path
        if self.option("sample") in ["column"]:
            path1 = self.work_dir + '/sample.xls'
            if os.path.exists(path1):
                os.remove(path1)
            os.link(sample, self.work_dir + '/sample.xls')
        else:
            with open(sample, 'r') as f:
                sample_list = [i.rstrip().split('\t') for i in f.readlines()]
                new_sample_list = map(lambda *a: '\t'.join(a) + '\n', *sample_list)
                with open(self.work_dir + '/sample.xls', 'w') as fw1:
                    fw1.writelines(new_sample_list)

    def get_sub_group(self):
        group_table = self.option("mann_group").path
        input_table = self.work_dir + '/sample.xls'
        with open(input_table, 'r') as s:
            s1 = s.readlines()[0]
            s_list = s1.rstrip().split('\t')
        if self.option('compare_type') != "multi":
            g_list = [self.option("group_name1"), self.option("group_name2")]
        else:
            if self.option("group_name").find(";") >= 0:
                g_list = self.option("group_name").split(';')
            else:
                g_list = self.option("group_name").split(',')
        new_list = list(set(g_list))
        t_list = ["name", "group"]
        df = pd.DataFrame()
        df.insert(0, 0, t_list)
        n = 1
        group_table_list = []
        with open(group_table, 'r') as g:
            for i in g.readlines():
                list1 = i.rstrip().split('\t')
                group_table_list.append(list1[1])
                if list1[1] in new_list and list1[0] in s_list:
                    n += 1
                    df.insert((n - 1), (n - 1), list1)
                else:
                    pass
            new_df = df.T
            new_df.to_csv(self.work_dir + '/sub_group_table.xls', sep='\t', header=None, index=0)
            new_df1 = new_df.drop(0, axis=0, inplace=False)
            new_df1.to_csv(self.work_dir + '/mann_sub_group_table.xls', sep='\t', header=None, index=0)
        for l in g_list:
            if l not in group_table_list:
                raise OptionError('选择比较的分组不在分组文件中，请检查！')

    def run_mann(self):
        two_group_test(self.work_dir + "/sample.xls", self.work_dir + "/mann_sub_group_table.xls",
                       self.work_dir + '/mann_result.xls', self.work_dir + '/mann_boxfile.xls', "mann",
                       str(1 - self.option('mann_ci')), self.option('mann_type'),
                       self.option('correction'),self.option('method'))
        cmd = self.r_path + " run_mann_test.r"
        self.logger.info("开始运行mann检验")
        command = self.add_command("mann_cmd", cmd, ignore_error=True).run()
        self.wait(command)
        def success():
            self.logger.info("生成单物种柱状图的数据")
            group_bar_toollab(self.work_dir + '/sample.xls', self.work_dir + "/mann_sub_group_table.xls",
                      self.work_dir + '/tmp_group_bar.xls', 'mann',self.option('method'))
            cmd1 = self.r_path + " run_mann_bar.r"
            bar_cmd = self.add_command("bar_cmd", cmd1, ignore_error=True).run()
            self.wait(bar_cmd)
            self.logger.info("开始运行计算置信区间")
            bootstrap(self.work_dir + '/tmp_group_bar.xls', self.work_dir + "/sub_group_table.xls",
                      self.option('coverage'),self.option('method'))
            def child_success():
                self.logger.info("mann_test运行完成")
            def child_fail():
                self.set_error("bar_cmd运行出错!")
            check_command(self, bar_cmd, [0], [-9], child_success, child_fail)
        def fail():
            self.set_error("mann_cmd运行出错!")
        check_command(self, command, [0], [-9], success, fail)

    def run_kru(self):
        if self.option('method') == "T":
            method = "true"
        else:
            method = "false"
        mul_group_test_toollab(self.work_dir + '/sample.xls', self.work_dir + '/kru_H_result.xls',
                       self.work_dir + '/kru_H_boxfile.xls', self.work_dir + "/sub_group_table.xls",
                       "kru_H", self.option('correction'),method,method)
        cmd = self.r_path + " run_kru_H_test.r"
        self.logger.info("开始运行mann检验")
        command = self.add_command("kru_h_cmd", cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("kru_cmd运行完成，开始运行post-hoc检验")
            in_file = self.work_dir + '/kru_H_result.xls'
            if self.option("kru_methor") in ["dunn", "nemenyi", "conover-iman", "steel_dwass"]:
                in_file = self.work_dir + '/sample.xls'
            self.posthoc(self.option("kru_methor"), in_file,
                         self.work_dir + "/sub_group_table.xls", self.option("coverage"), './kru_H')
            self.logger.info("生成单物种柱状图的数据")
            group_bar_toollab(self.work_dir + '/sample.xls', self.work_dir + "/sub_group_table.xls",
                      self.work_dir + '/tmp_group_bar.xls', 'kru_H',self.option('method'))
            cmd1 = self.r_path + " run_kru_H_bar.r"
            bar_cmd = self.add_command("bar_cmd", cmd1).run()
            self.wait(bar_cmd)
            self.logger.info("开始运行计算置信区间")
            bootstrap(self.work_dir + '/tmp_group_bar.xls', self.work_dir + "/sub_group_table.xls",
                      self.option('coverage'),self.option('method'))
            if bar_cmd.return_code == 0:
                self.logger.info("kru_H_test运行完成")
            else:
                self.set_error("bar_cmd运行出错!")
        else:
            self.set_error("kru_cmd运行出错!")

    def posthoc(self, methor, statfile, groupfile, coverage, outfile):
        if methor == 'tukeykramer':
            tukeykramer(statfile, groupfile, coverage, outfile)
        elif methor == 'gameshowell':
            gameshowell(statfile, groupfile, coverage, outfile)
        elif methor == 'welchuncorrected':
            welchuncorrected(statfile, groupfile, coverage, outfile)
        elif methor == 'scheffe':
            scheffe(statfile, groupfile, coverage, outfile)
        else:
            self.logger.info("进行kru独有的 posthoc 计算{}".format(methor))
            r_path = "/bioinfo/metaGenomic/metaphlan3/bin/Rscript"
            r_posthoc = self.config.PACKAGE_DIR + '/statistical/posthoc.r'
            outfile = "{}_{}.xls".format(outfile, methor)
            cmd = "{} {} {} {} {} {}".format(r_path, r_posthoc, statfile, groupfile, methor, outfile)
            command = self.add_command('r_posthoc', cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("kru posthoc 计算完成")
            else:
                self.set_error("kru posthoc 计算失败")

    def get_box(self):
        box_file = self.work_dir + '/box_result.xls'
        get_box_value_row(self.option('mann_input').path, self.work_dir + "/sub_group_table.xls", box_file,self.option('method'))

    def link_table(self, infile_path, outfile_path):
        with open(outfile_path, 'w'):
            path1 = outfile_path
            if os.path.exists(path1):
                os.remove(path1)
            os.link(infile_path, outfile_path)

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        """
        if self.option('compare_type') != "multi":
            re_df = pd.read_table('mann_result.xls')
            re_df.drop(['statistic'], axis=1, inplace=True)
            col_list_re = re_df.columns.tolist()
            re_df.rename(columns={col_list_re[0]: "Name"}, inplace=True)
            ci_df = pd.read_table('bootstrap_CI.xls')
            col_list_ci = ci_df.columns.tolist()
            ci_df.rename(columns={col_list_ci[0]: "Name"}, inplace=True)
            new_re = pd.merge(re_df, ci_df, on="Name")
            new_re.rename(columns={'lowerCI': "Lower ci", "upperCI": "Upper ci"}, inplace=True)
            #new_re.drop(['effectsize'], axis=1, inplace=True)
            new_re_df = pd.DataFrame()
            #effectsize_list = []
            #for j in range(len(col_list_re)):
            #    if col_list_re[j].find('mean') != -1:
            #        effectsize_list.append(col_list_re[j])
            #        new_re_df[col_list_re[j]] = re_df[col_list_re[j]]
            new_re_df['Effectsize'] = new_re["effectsize"]
            new_re.drop(['effectsize'], axis=1, inplace=True)
            new_re['Effectsize'] = new_re_df['Effectsize']
            new_re.to_csv(self.work_dir + '/compare_result0.xls', sep='\t')
            self.link_table(self.work_dir + '/compare_result0.xls', self.output_dir + '/compare_result.xls')
            self.link_table(self.work_dir + '/mann_boxfile.xls', self.output_dir + '/compare_boxfile.xls')
            self.link_table(self.work_dir + '/bootstrap_CI.xls', self.output_dir + '/compare_CI.xls')
            self.link_table(self.work_dir + '/tmp_group_bar.xls', self.output_dir + '/tmp_group_bar.xls')
        else:
            if ";" in self.option("group_name"):
                g_list = self.option("group_name").split(';')
            elif "," in self.option("group_name"):
                g_list = self.option("group_name").split(',')
            # a = g_list[0]
            # b = g_list[1]
            # c = g_list[2]
            m = self.option("kru_methor")
            re_df = pd.read_table(self.work_dir + '/kru_H_result.xls')
            re_df.drop(['statistic'], axis=1, inplace=True)
            re_df.rename(columns={"qvalue": "corrected_pvalue"}, inplace=True)
            col_list = re_df.columns.tolist()
            in_list = list(re_df.index)
            re_df.index = range(len(in_list))
            re_df["Name"] = in_list
            order = col_list
            order.insert(0, 'Name')
            re_df = re_df[order]
            new_re_df = pd.DataFrame()
            for j in range(len(col_list)):
                if col_list[j].find('mean') != -1:
                    new_re_df[col_list[j]] = re_df[col_list[j]]
            new_re_df['sum'] = new_re_df.apply(lambda x: x.sum(), axis=1)
            re_df['Abundance'] = new_re_df['sum']
            re_df.to_csv(self.work_dir + '/compare_result0.xls', sep='\t')
            self.link_table(self.work_dir + '/compare_result0.xls', self.output_dir + '/compare_result.xls')
            self.link_table(self.work_dir + '/kru_H_boxfile.xls', self.output_dir + '/compare_boxfile.xls')
            self.link_table(self.work_dir + '/bootstrap_CI.xls', self.output_dir + '/compare_CI.xls')
            self.link_table(self.work_dir + '/tmp_group_bar.xls', self.output_dir + '/tmp_group_bar.xls')
            
            kru_H_files = glob.glob(os.path.join(self.work_dir + '/', 'kru_H' + '_' + m + '_*.xls'))
            print (kru_H_files)
            if len(kru_H_files) > 0:
                for k_file in kru_H_files:
                    k_file_name = os.path.basename(k_file)
                    self.link_table(k_file, self.output_dir + '/' + k_file_name)
            else:
                kru_H_files = glob.glob(os.path.join(self.work_dir + '/', 'kru_H' + '_' + m + '.xls'))
                print (kru_H_files)
                if len(kru_H_files) > 0:
                    for k_file in kru_H_files:
                        k_file_name = os.path.basename(k_file)
                        self.link_table(k_file, self.output_dir + '/' + k_file_name)
            # if os.path.exists(self.work_dir + '/kru_H' + '_' + m + '_' + a + '-' + b + '.xls'):
            #     self.link_table(self.work_dir + '/kru_H' + '_' + m + '_' + a + '-' + b + '.xls',
            #                     self.output_dir + '/mul1_kru_H.xls')
            # else:
            #     self.link_table(self.work_dir + '/kru_H' + '_' + m + '_' + b + '-' + a + '.xls',
            #                     self.output_dir + '/mul1_kru_H.xls')
            # if os.path.exists(self.work_dir + '/kru_H' + '_' + m + '_' + b + '-' + c + '.xls'):
            #     self.link_table(self.work_dir + '/kru_H' + '_' + m + '_' + b + '-' + c + '.xls',
            #                     self.output_dir + '/mul2_kru_H.xls')
            # else:
            #     self.link_table(self.work_dir + '/kru_H' + '_' + m + '_' + c + '-' + b + '.xls',
            #                     self.output_dir + '/mul2_kru_H.xls')
            # if os.path.exists(self.work_dir + '/kru_H' + '_' + m + '_' + a + '-' + c + '.xls'):
            #     self.link_table(self.work_dir + '/kru_H' + '_' + m + '_' + a + '-' + c + '.xls',
            #                     self.output_dir + '/mul3_kru_H.xls')
            # else:
            #     self.link_table(self.work_dir + '/kru_H' + '_' + m + '_' + c + '-' + a + '.xls',
            #                     self.output_dir + '/mul3_kru_H.xls')
        self.link_table(self.work_dir + '/box_result.xls', self.output_dir + '/box_result_group.xls')
