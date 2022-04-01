# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.taxon.mask_taxon import mask_taxon
from mbio.packages.metabolome.scripts.metastat import *
# from mbio.packages.statistical.metastat import *
from mbio.packages.statistical.twogroup_CI import *
from mbio.packages.statistical.twosample_CI import *
# from mbio.packages.statistical.mul_posthoc import *
from mbio.packages.metabolome.scripts.mul_posthoc import *
import os
import re
import pandas as pd


class MetastatAgent(Agent):
    """
    statistical metastat+ 调用metastat.py 包进行差异性分析
    version v1.0
    author: qiuping
    last_modify: 2015.11.30
    """
    def __init__(self, parent):
        super(MetastatAgent, self).__init__(parent)
        options = [
            {"name": "chi_input", "type": "infile", "format": "meta.otu.otu_table"},  # 卡方检验的输入文件
            {"name": "chi_sample1", "type": "string"},  # 卡方检验的输入样品名称
            {"name": "chi_sample2", "type": "string"},  # 卡方检验的输入样品名称
            {"name": "chi_correction", "type": "string", "default": "none"},  # 卡方检验的多重检验校正
            {"name": "fisher_input", "type": "infile", "format": "meta.otu.otu_table"},  # 费舍尔检验的输入文件
            {"name": "fisher_ci", "type": "float", "default": 0.05},  # 费舍尔检验的显著性水平
            {"name": "fisher_sample1", "type": "string"},  # 费舍尔检验的输入样品名称1
            {"name": "fisher_sample2", "type": "string"},  # 费舍尔检验的输入样品名称2
            {"name": "fisher_correction", "type": "string", "default": "none"},  # 费舍尔检验的多重检验校正
            {"name": "fisher_type", "type": "string", "default": "two.side"},  # 费舍尔检验的选择单尾或双尾检验
            {"name": "kru_H_input", "type": "infile", "format": "meta.otu.otu_table"},  # kruskal_wallis_H_test的输入文件
            {"name": "kru_H_group", "type": "infile", "format": "meta.otu.group_table,toolapps.group_table"},  # kruskal_wallis_H_test输入分组
            # {"name": "kru_H_type", "type": "string", "default": "two.side"},  #kruskal_wallis_H_test选择单双尾检验
            {"name": "kru_H_correction", "type": "string", "default": "none"},  # kruskal_wallis_H_test的多重检验校正
            {"name": "mann_input", "type": "infile", "format": "meta.otu.otu_table"},  # 秩和检验的输入文件
            {"name": "mann_ci", "type": "float", "default": 0.05},  # 秩和检验的显著性水平
            {"name": "mann_group", "type": "infile", "format": "meta.otu.group_table,toolapps.group_table"},  # 秩和检验的输入分组文件
            {"name": "mann_correction", "type": "string", "default": "none"},  # 秩和检验的多重检验校正
            {"name": "mann_type", "type": "string", "default": "two.side"},  # 秩和检验的选择单尾或双尾检验
            {"name": "signal_input", "type": "infile", "format": "meta.otu.otu_table"},  # 符号秩和检验的输入文件
            {"name": "signal_ci", "type": "float", "default": 0.05},  # 符号秩和检验的显著性水平
            {"name": "signal_group", "type": "infile", "format": "meta.otu.group_table,toolapps.group_table"},  # 符号秩和检验的输入分组文件
            {"name": "signal_correction", "type": "string", "default": "none"},  # 符号秩和检验的多重检验校正
            {"name": "signal_type", "type": "string", "default": "two.side"},  # 符号秩和检验的选择单尾或双尾检验
            {"name": "student_input", "type": "infile", "format": "meta.otu.otu_table"},  # T检验的输入文件
            {"name": "student_ci", "type": "float", "default": 0.05},  # T检验的显著性水平
            {"name": "student_group", "type": "infile", "format": "meta.otu.group_table,toolapps.group_table"},  # T检验的输入分组文件
            {"name": "student_correction", "type": "string", "default": "none"},  # T检验的多重检验校正
            {"name": "student_type", "type": "string", "default": "two.side"},  # T检验的选择单尾或双尾检验
            {"name": "welch_input", "type": "infile", "format": "meta.otu.otu_table"},  # welch_T检验的输入文件
            {"name": "welch_ci", "type": "float", "default": 0.05},  # welch_T检验的显著性水平
            {"name": "welch_group", "type": "infile", "format": "meta.otu.group_table,toolapps.group_table"},  # welch_T检验的输入分组文件
            {"name": "welch_correction", "type": "string", "default": "none"},  # welch_T检验的多重检验校正
            {"name": "welch_type", "type": "string", "default": "two.side"},  # welch_T检验的选择单尾或双尾检验
            {"name": "anova_input", "type": "infile", "format": "meta.otu.otu_table"},  # anova分析的输入文件
            {"name": "anova_group", "type": "infile", "format": "meta.otu.group_table,toolapps.group_table"},  # anova分析的输入分组文件
            {"name": "anova_correction", "type": "string", "default": "none"},  # anova分析的多重检验校正
            {"name": "test", "type": "string"},   # 选择统计学检验分析方法
            {"name": "student_gname", "type": "string"},  # student检验分组方案选择
            {"name": "welch_gname", "type": "string"},  # welch检验分组方案选择
            {"name": "mann_gname", "type": "string"},  # wilcox秩和检验分组方案选择
            {"name": "signal_gname", "type": "string"},  # 符号秩和检验分组方案选择
            {"name": "kru_H_gname", "type": "string"},  # kru检验分组方案选择
            {"name": "anova_gname", "type": "string"},  # 单因素方差分析分组方案选择
            {"name": "kru_H_coverage", "type": "float", "default": 0.95},  # 计算置信区间所选择的置信度
            {"name": "anova_coverage", "type": "float", "default": 0.95},
            {"name": "student_coverage", "type": "float", "default": 0.95},
            {"name": "welch_coverage", "type": "float", "default": 0.95},
            {"name": "mann_coverage", "type": "float", "default": 0.95},
            {"name": "signal_coverage", "type": "float", "default": 0.95},
            {"name": "chi_coverage", "type": "float", "default": 0.95},
            {"name": "fisher_coverage", "type": "float", "default": 0.95},
            {"name": "kru_H_methor", "type": "string", "default": 'tukeykramer'},  # post-hoc检验的方法
            {"name": "anova_methor", "type": "string", "default": 'tukeykramer'},  # post-hoc检验的方法
            {"name": "chi_methor", "type": "string", "default": 'DiffBetweenPropAsymptotic'},  # 两样本计算置信区间的方法
            {"name": "fisher_methor", "type": "string", "default": 'DiffBetweenPropAsymptotic'},  # 两样本计算置信区间的方法
            {"name": "est_group", "type": "infile", "format": "meta.alpha_diversity.group_file_dir"},
            {"name": "est_input", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "est_test_method", "type": "string", "default": "student"}  # add by qindanhua 20170105
        ]
        self.add_option(options)
        self.step.add_steps("stat_test")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.stat_test.start()
        self.step.update()

    def stepfinish(self):
        self.step.stat_test.finish()
        self.step.update()

    def check_sample(self, group_sample, otu_sample):
        for name in group_sample:
            if name not in otu_sample:
                raise OptionError('分组样本不在otu表所拥有的样本内，请检查分组方案', code="32502901")

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
                raise OptionError('一个分组类别下面至少有两个样本，请重新分组', code="32502902")
            if len(sample) != len(set(sample)):
                raise OptionError('不同的分组类别下有相同的样本，此分析不支持该分组方案，请重新分组', code="32502903")

    def check_options(self):
        """
        检查参数设置
        :return:
        """
        if not self.option('test'):
            self.logger.info(self.option('test'))
            raise OptionError("必须设置输入的检验名称", code="32502904")
        for i in self.option('test').split(','):
            if i not in ["chi", "fisher", "kru_H", "mann", "signal", "anova", "student", "welch", "estimator"]:
                raise OptionError("所输入的检验名称不对", code="32502905")
            elif i == "chi":
                if not self.option("chi_input").is_set:
                    raise OptionError('必须设置卡方检验输入的otutable文件', code="32502906")
                if not self.option("chi_sample1") and not self.option("chi_sample2"):
                    raise OptionError('必须设置卡方检验要比较的样品名', code="32502907")
                if self.option("chi_correction") not in ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
                                                         "none"]:
                    raise OptionError('卡方检验的多重检验校正的方法不被支持', code="32502908")
                if self.option("chi_coverage") not in [0.90, 0.95, 0.98, 0.99, 0.999]:
                    raise OptionError('卡方检验的置信区间的置信度不在范围值内', code="32502909")
                if self.option("chi_methor") not in ["DiffBetweenPropAsymptoticCC", "DiffBetweenPropAsymptotic",
                                                     "NewcombeWilson"]:
                    raise OptionError('卡方检验的计算置信区间的方法不在范围值内', code="32502910")
                otu_sample = self.option('chi_input').get_sample_info()
                if self.option('chi_sample1') not in otu_sample or self.option('chi_sample2') not in otu_sample:
                    raise OptionError('输入的样本不在otu表所拥有的样本内，请检查卡方检验样本名', code="32502911")
            elif i == "fisher":
                if not self.option("fisher_input").is_set:
                    raise OptionError('必须设置费舍尔检验输入的otutable文件', code="32502912")
                if not self.option("fisher_sample1") and not self.option("fisher_sample2"):
                    raise OptionError('必须设置费舍尔检验要比较的样品名', code="32502913")
                if self.option("fisher_correction") not in ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                                                            "fdr", "none"]:
                    raise OptionError('该多重检验校正的方法不被支持', code="32502914")
                if self.option("fisher_ci") <= 0 or self.option("fisher_ci") >= 1:
                    raise OptionError('所输入的显著水平不在范围值内', code="32502915")
                if self.option("fisher_type") not in ["two.side", "greater", "less"]:
                    raise OptionError('所输入的类型不在范围值内', code="32502916")
                if self.option("fisher_coverage") not in [0.90, 0.95, 0.98, 0.99, 0.999]:
                    raise OptionError('费舍尔检验的置信区间的置信度不在范围值内', code="32502917")
                if self.option("fisher_methor") not in ["DiffBetweenPropAsymptoticCC", "DiffBetweenPropAsymptotic",
                                                        "NewcombeWilson"]:
                    raise OptionError('费舍尔检验的计算置信区间的方法不在范围值内', code="32502918")
                otu_sample = self.option('fisher_input').get_sample_info()
                if self.option('fisher_sample1') not in otu_sample or self.option('fisher_sample2') not in otu_sample:
                    raise OptionError('输入的样本不在otu表所拥有的样本内，请检查费舍尔检验样本名', code="32502919")
            elif i == "kru_H":
                if not self.option("kru_H_input").is_set:
                    raise OptionError('必须设置Kruskal-Wallis秩和检验输入的otutable文件', code="32502920")
                if not self.option("kru_H_group").is_set:
                    raise OptionError('必须设置Kruskal-Wallis秩和检验输入的分组文件', code="32502921")
                if not self.option("kru_H_gname"):
                    raise OptionError("Kruskal-Wallis秩和检验分组方案选择参数为必须参数，请设置", code="32502922")
                if len(self.option("kru_H_gname").split(',')) != 1:
                    raise OptionError("组间差异的分组方案只能为1个", code="32502923")
                if self.option("kru_H_correction") not in ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                                                           "fdr", "none"]:
                    raise OptionError('该多重检验校正的方法不被支持', code="32502924")
                # if self.option("kru_H_type") not in ["two.side", "greater", "less"]:
                #     raise OptionError('所输入的类型不在范围值内')
                gnum = self.option('kru_H_group').group_num(self.option('kru_H_gname'))
                if gnum < 3:
                    raise OptionError("Kruskal-Wallis秩和检验的分组方案的分组类别必须大于等于3", code="32502925")
                if self.option("kru_H_coverage") not in [0.90, 0.95, 0.98, 0.99, 0.999]:
                    raise OptionError('Kruskal-Wallis秩和检验的posthoc的置信度不在范围值内', code="32502926")
                if self.option("kru_H_methor") not in ["scheffe", "welchuncorrected", "tukeykramer", "gameshowell"]:
                    raise OptionError('Kruskal-Wallis秩和检验的posthoc检验方法不在范围值内', code="32502927")
                otu_sample = self.option('kru_H_input').get_sample_info()
                self.logger.info(otu_sample)
                group_sample, header, is_empty = self.option('kru_H_group').get_file_info()
                self.check_sample(group_sample, otu_sample)
                self.check_group(self.option('kru_H_group').prop['path'])
            elif i == "anova":
                if not self.option("anova_input").is_set:
                    raise OptionError('必须设置Kruskal-Wallis秩和检验输入的otutable文件', code="32502928")
                if not self.option("anova_group").is_set:
                    raise OptionError('必须设置Kruskal-Wallis秩和检验输入的分组文件', code="32502929")
                if not self.option("anova_gname"):
                    raise OptionError("单因素方差分析分组方案选择参数为必须参数，请设置", code="32502930")
                if len(self.option("anova_gname").split(',')) != 1:
                    raise OptionError("组间差异的分组方案只能为1个", code="32502931")
                if self.option("anova_correction") not in ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                                                           "fdr", "none"]:
                    raise OptionError('该多重检验校正的方法不被支持', code="32502932")
                gnum = self.option('anova_group').group_num(self.option('anova_gname'))
                if gnum < 3:
                    raise OptionError("单因素方差分析的分组方案的分组类别必须大于等于3", code="32502933")
                if self.option("anova_coverage") not in [0.90, 0.95, 0.98, 0.99, 0.999]:
                    raise OptionError('单因素方差分析的posthoc的置信度不在范围值内', code="32502934")
                if self.option("anova_methor") not in ["scheffe", "welchuncorrected", "tukeykramer", "gameshowell"]:
                    raise OptionError('单因素方差分析的posthoc检验方法不在范围值内', code="32502935")
                otu_sample = self.option('anova_input').get_sample_info()
                self.logger.info(otu_sample)
                group_sample, header, is_empty = self.option('anova_group').get_file_info()
                self.check_sample(group_sample, otu_sample)
                self.check_group(self.option('anova_group').prop['path'])
            elif i == "mann":
                if not self.option("mann_input").is_set:
                    raise OptionError('必须设置wilcox秩和检验输入的otutable文件', code="32502936")
                if not self.option("mann_group").is_set:
                    raise OptionError('必须设置wilcox秩和检验输入的分组文件', code="32502937")
                if self.option("mann_correction") not in ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
                                                          "none"]:
                    raise OptionError('该多重检验校正的方法不被支持', code="32502938")
                if self.option("mann_ci") <= 0 or self.option("mann_ci") >= 1:
                    raise OptionError('所输入的显著水平不在范围值内', code="32502939")
                if self.option("mann_type") not in ["two.side", "greater", "less"]:
                    raise OptionError('所输入的类型不在范围值内', code="32502940")
                if not self.option("mann_gname"):
                    raise OptionError("wilcox秩和检验分组方案选择参数为必须参数，请设置", code="32502941")
                if len(self.option("mann_gname").split(',')) != 1:
                    raise OptionError("组间差异的分组方案只能为1个", code="32502942")
                gnum = self.option('mann_group').group_num(self.option('mann_gname'))
                if gnum != 2:
                    raise OptionError("wilcox秩和检验的分组方案的分组类别必须等于2", code="32502943")
                if self.option("mann_coverage") not in [0.90, 0.95, 0.98, 0.99, 0.999]:
                    raise OptionError('wilcox秩和检验的置信区间的置信度不在范围值内', code="32502944")
                otu_sample = self.option('mann_input').get_sample_info()
                self.logger.info(otu_sample)
                group_sample, header, is_empty = self.option('mann_group').get_file_info()
                self.check_sample(group_sample, otu_sample)
                self.check_group(self.option('mann_group').prop['path'])
            elif i == "signal":
                if not self.option("signal_input").is_set:
                    raise OptionError('必须设置wilcox秩和检验输入的otutable文件', code="32502945")
                if not self.option("signal_group").is_set:
                    raise OptionError('必须设置wilcox秩和检验输入的分组文件', code="32502946")
                if self.option("signal_correction") not in ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
                                                          "none"]:
                    raise OptionError('该多重检验校正的方法不被支持', code="32502947")
                if self.option("signal_ci") <= 0 or self.option("signal_ci") >= 1:
                    raise OptionError('所输入的显著水平不在范围值内', code="32502948")
                if self.option("signal_type") not in ["two.side", "greater", "less"]:
                    raise OptionError('所输入的类型不在范围值内', code="32502949")
                if not self.option("signal_gname"):
                    raise OptionError("wilcox秩和检验分组方案选择参数为必须参数，请设置", code="32502950")
                if len(self.option("signal_gname").split(',')) != 1:
                    raise OptionError("组间差异的分组方案只能为1个", code="32502951")
                gnum = self.option('signal_group').group_num(self.option('signal_gname'))
                if gnum != 2:
                    raise OptionError("wilcox秩和检验的分组方案的分组类别必须等于2", code="32502952")
                if self.option("signal_coverage") not in [0.90, 0.95, 0.98, 0.99, 0.999]:
                    raise OptionError('wilcox秩和检验的置信区间的置信度不在范围值内', code="32502953")
                otu_sample = self.option('signal_input').get_sample_info()
                self.logger.info(otu_sample)
                group_sample, header, is_empty = self.option('signal_group').get_file_info()
                self.check_sample(group_sample, otu_sample)
                self.check_group(self.option('signal_group').prop['path'])
            elif i == "student":
                if not self.option("student_input").is_set:
                    raise OptionError('必须设置student_T检验输入的otutable文件', code="32502954")
                if not self.option("student_group").is_set:
                    raise OptionError('必须设置student_T检验输入的分组文件', code="32502955")
                if self.option("student_correction") not in ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                                                             "fdr", "none"]:
                    raise OptionError('该多重检验校正的方法不被支持', code="32502956")
                if self.option("student_ci") <= 0 or self.option("student_ci") >= 1:
                    raise OptionError('所输入的显著水平不在范围值内', code="32502957")
                if self.option("student_type") not in ["two.side", "greater", "less"]:
                    raise OptionError('所输入的类型不在范围值内', code="32502958")
                if not self.option("student_gname"):
                    raise OptionError("student_T检验分组方案选择参数为必须参数，请设置", code="32502959")
                if len(self.option("student_gname").split(',')) != 1:
                    raise OptionError("组间差异的分组方案只能为1个", code="32502960")
                gnum = self.option('student_group').group_num(self.option('student_gname'))
                if gnum != 2:
                    raise OptionError("student_T检验的分组方案的分组类别必须等于2", code="32502961")
                if self.option("student_coverage") not in [0.90, 0.95, 0.98, 0.99, 0.999]:
                    raise OptionError('student_T检验的置信区间的置信度不在范围值内', code="32502962")
                otu_sample = self.option('student_input').get_sample_info()
                self.logger.info(otu_sample)
                group_sample, header, is_empty = self.option('student_group').get_file_info()
                self.check_sample(group_sample, otu_sample)
                self.check_group(self.option('student_group').prop['path'])
            elif i == "welch":
                if not self.option("welch_input").is_set:
                    raise OptionError('必须设置welch_T检验输入的otutable文件', code="32502963")
                if not self.option("welch_group").is_set:
                    raise OptionError('必须设置welch_T检验输入的分组文件', code="32502964")
                if self.option("welch_correction") not in ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                                                           "fdr", "none"]:
                    raise OptionError('该多重检验校正的方法不被支持', code="32502965")
                if self.option("welch_ci") <= 0 or self.option("welch_ci") >= 1:
                    raise OptionError('所输入的显著水平不在范围值内', code="32502966")
                if self.option("welch_type") not in ["two.side", "greater", "less"]:
                    raise OptionError('所输入的类型不在范围值内', code="32502967")
                if not self.option("welch_gname"):
                    raise OptionError("welch_T检验分组方案选择参数为必须参数，请设置", code="32502968")
                if len(self.option("welch_gname").split(',')) != 1:
                    raise OptionError("组间差异的分组方案只能为1个", code="32502969")
                gnum = self.option('welch_group').group_num(self.option('welch_gname'))
                if gnum != 2:
                    raise OptionError("welch_T检验的分组方案的分组类别必须等于2", code="32502970")
                if self.option("welch_coverage") not in [0.90, 0.95, 0.98, 0.99, 0.999]:
                    raise OptionError('welch_T检验的置信区间的置信度不在范围值内', code="32502971")
                otu_sample = self.option('welch_input').get_sample_info()
                self.logger.info(otu_sample)
                group_sample, header, is_empty = self.option('welch_group').get_file_info()
                self.check_sample(group_sample, otu_sample)
                self.check_group(self.option('welch_group').prop['path'])
            elif i == "estimator":
                if not self.option("est_input").is_set:
                    raise OptionError('必须设置est_T检验输入的多样性指数文件', code="32502972")
                if not self.option("est_group").is_set:
                    raise OptionError('必须设置est_T检验输入的分组文件夹', code="32502973")
                # otu_sample = self.option('est_input').get_sample_info()
                # self.logger.info(otu_sample)
                # group_sample = self.option('est_group').get_file_info()
                # self.check_sample(group_sample, otu_sample)
        return True

    def set_resource(self):
        """
        设置所需资源
        :return:
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            [r".*_result\.xls", "xls", "物种组间差异显著性比较结果表，包括均值，标准差，p值"],
            [r".*_CI\.xls", "xls", "组间差异显著性比较两组，两样本比较的置信区间值以及效果量"],
            [r".*(-).*", "xls", "组间差异显著性比较多组比较的posthoc检验比较的结果，包含置信区间，效果量，p值"],
            [r".*_boxfile\.xls", "xls", "组间差异显著性比较用于画箱线图的数据，包含四分位值"],
            [r"^est_result", "xls", "多样性指数T检验结果分析表，包括均值，标准差，p值，q值"]
            ])
        super(MetastatAgent, self).end()


class MetastatTool(Tool):
    def __init__(self, config):
        super(MetastatTool, self).__init__(config)
        self._version = "v1.0.1"
        self.package_path = 'packages/metastat.py'
        self.r_path = '/program/R-3.3.1/bin/Rscript'
        self.name_to_name = {}

    def remove_zero_and_mask_taxon(self, profile, group):
        if type(group) == str:
            samp_list = open(group).readlines()[1:]
            for i in range(len(samp_list)):
                samp_list[i] = samp_list[i].split('\t')[0]
        elif type(group) == list:
            samp_list = group
        else:
            raise OptionError('函数中group类型问题，type is %s' , variables=( type(group)), code="32502974")
        profile.sub_otu_sample(samp_list, self.work_dir + '/no_zero.xls')
        profile_path = self.work_dir + '/tmp_mask_otu.xls'
        self.name_to_name = mask_taxon(self.work_dir + '/no_zero.xls', profile_path)
        if len(self.name_to_name) < 2:
            self.set_error("所含非零物种/功能过少，更改其他分类水平分析。非零物种/功能数量为：%s" , variables=( len(self.name_to_name)), code="32502901")

    def dashrepl(self, matchobj):
        """
        add func by guhaidong 20171031
        """
        return self.name_to_name[matchobj.groups()[0]]

    def add_taxon(self, old_result, taxon_result):
        """
        add func by guhaidong 20171031
        description: 将旧注释的名称，根据词典替换成新注释名称
        """
        with open(old_result, "r") as f, open(taxon_result, "w") as w:
            # w.write(old_result)
            for i in f.readlines():
                #line = i.strip()
                new_line = re.sub(r"(name\d+)", self.dashrepl, i)
                w.write(new_line)

    def run(self):
        """
        运行
        :return:
        """
        super(MetastatTool, self).run()
        self.run_test()
        self.set_output()
        self.end()

    def run_test(self):
        """
        运行metastat.py
        :return:
        """
        for test in self.option('test').split(','):
            # self.logger.info(t)
            if test == "chi":
                self.run_chi()
            elif test == "fisher":
                self.run_fisher()
            elif test == "student":
                self.run_student()
            elif test == "welch":
                self.run_welch()
            elif test == "mann":
                self.run_mann()
            elif test == "signal":
                self.run_signal()
            elif test == "kru_H":
                self.run_kru()
            elif test == "anova":
                self.run_anova()
            elif test == "estimator":
                self.run_est()

    def run_est(self):    # modify by qindanhua 20170105 add option test method
        gfilelist = os.listdir(self.option("est_group").prop['path'])
        self.logger.info(gfilelist)
        i = 1
        for group in gfilelist:
            est_ttest(self.option('est_input').prop['path'], self.work_dir + '/est_result%s.xls' % i,
                      os.path.join(self.option("est_group").prop['path'], group), self.option("est_test_method"))
            cmd = self.r_path + " run_est_ttest.r"
            self.logger.info("开始运行est_T检验")
            command = self.add_command("est_cmd{}".format(i), cmd).run()
            i += 1
            self.wait(command)
            if command.return_code == 0 or command.return_code is None:
                self.logger.info("est_ttest运行完成")
            else:
                self.set_error("est_ttest运行出错!", code="32502902")

    def run_chi(self):
        self.remove_zero_and_mask_taxon(self.option('chi_input'), [self.option('chi_sample1'), self.option('chi_sample2')])
        # self.name_to_name = mask_taxon(self.option('chi_input').prop['path'], self.work_dir + "/tmp_mask_otu.xls")
        two_sample_test(self.work_dir + "/tmp_mask_otu.xls", self.work_dir + '/chi_result.xls', "chi",
                        self.option('chi_sample1'), self.option('chi_sample2'), self.option('chi_correction'))
        # two_sample_test(self.option('chi_input').prop['path'], self.work_dir + '/chi_result.xls', "chi",
        #                 self.option('chi_sample1'), self.option('chi_sample2'), self.option('chi_correction'))
        cmd = self.r_path + " run_chi_test.r"
        self.logger.info("开始运行卡方检验")
        command = self.add_command("chi_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("chi_cmd运行完成，开始运行计算置信区间")
            # self.twosample_ci(self.option("chi_methor"), self.option("chi_input").prop['path'],
            self.twosample_ci(self.option("chi_methor"), self.work_dir + "/tmp_mask_otu.xls",
                              self.work_dir + '/chi_result.xls', self.option('chi_sample1'), self.option('chi_sample2'),
                              self.option('chi_coverage'), self.work_dir + '/chi_CI.xls')
        else:
            self.set_error("chi_cmd运行出错!", code="32502903")
        self.logger.info("chi_test运行完成")

    def run_fisher(self):
        self.remove_zero_and_mask_taxon(self.option('fisher_input'), [self.option('fisher_sample1'), self.option('fisher_sample2')])
        # self.name_to_name = mask_taxon(self.option("fisher_input").prop['path'], self.work_dir + "/tmp_mask_otu.xls")
        two_sample_test(self.work_dir + '/tmp_mask_otu.xls', self.work_dir + '/fisher_result.xls', "fisher",
                        self.option('fisher_sample1'), self.option('fisher_sample2'),
                        str(1 - self.option('fisher_ci')), self.option('fisher_type'),
                        self.option('fisher_correction'))
        # two_sample_test(self.option('fisher_input').prop['path'], self.work_dir + '/fisher_result.xls', "fisher",
        #                 self.option('fisher_sample1'), self.option('fisher_sample2'),
        #                 str(1 - self.option('fisher_ci')), self.option('fisher_type'),
        #                 self.option('fisher_correction'))
        cmd = self.r_path + " run_fisher_test.r"
        self.logger.info("开始运行fisher检验")
        command = self.add_command("fisher_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("fisher_cmd运行完成，开始运行计算置信区间")
            self.twosample_ci(self.option("fisher_methor"), self.work_dir + '/tmp_mask_otu.xls',
                              self.work_dir + '/fisher_result.xls', self.option('fisher_sample1'),
                              self.option('fisher_sample2'), self.option('fisher_coverage'),
                              self.work_dir + '/fisher_CI.xls')
            # self.twosample_ci(self.option("fisher_methor"), self.option("fisher_input").prop['path'],
            #                   self.work_dir + '/fisher_result.xls', self.option('fisher_sample1'),
            #                   self.option('fisher_sample2'), self.option('fisher_coverage'),
            #                   self.work_dir + '/fisher_CI.xls')
        else:
            self.set_error("fisher_cmd运行出错!", code="32502904")
        self.logger.info("fisher_test运行完成")

    def twosample_ci(self, methor, otufile, statfile, sample1, sample2, coverage, outfile):
        if methor == "DiffBetweenPropAsymptoticCC":
            DiffBetweenPropAsymptoticCC(otufile, statfile, sample1, sample2, coverage, outfile)
        if methor == "DiffBetweenPropAsymptotic":
            DiffBetweenPropAsymptotic(otufile, statfile, sample1, sample2, coverage, outfile)
        if methor == "NewcombeWilson":
            NewcombeWilson(otufile, statfile, sample1, sample2, coverage, outfile)

    def run_student(self):
        glist = [self.option('student_gname')]
        self.option('student_group').sub_group('./student_group', glist)
        self.remove_zero_and_mask_taxon(self.option('student_input'), './student_group')
        two_group_test(self.work_dir + "/tmp_mask_otu.xls", './student_group',
        # two_group_test(self.option('student_input').prop['path'], './student_group',
                       self.work_dir + '/student_result.xls', self.work_dir + '/student_boxfile.xls', "student",
                       str(1 - self.option('student_ci')), self.option('student_type'),
                       self.option('student_correction'))
        cmd = self.r_path + " run_student_test.r"
        self.logger.info("开始运行student_T检验")
        command = self.add_command("student_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("student_cmd运行完成，开始运行计算置信区间")
            student(self.work_dir + '/student_result.xls', './student_group', self.option('student_coverage'))
            self.logger.info("生成单物种柱状图的数据")
            group_bar(self.work_dir + '/tmp_mask_otu.xls', './student_group', self.work_dir + '/tmp_group_bar.xls', 'student')
            # group_bar(self.option('student_input').prop['path'], './student_group', self.work_dir + '/student_plot_group_bar.xls', "student")
            cmd1 = self.r_path + " run_student_bar.r"
            bar_cmd = self.add_command("bar_cmd", cmd1).run()
            self.wait(bar_cmd)
            self.add_taxon(self.work_dir + '/tmp_group_bar.xls', self.work_dir + '/student_plot_group_bar.xls')
            if bar_cmd.return_code == 0:
                self.logger.info("student_test运行完成")
            else:
                self.set_error("bar_cmd运行出错!", code="32502905")
        else:
            self.set_error("student_cmd运行出错!", code="32502906")

    def run_welch(self):
        # self.name_to_name = mask_taxon(self.option("welch_input").prop['path'], self.work_dir + "/tmp_mask_otu.xls")
        glist = [self.option('welch_gname')]
        self.option('welch_group').sub_group('./welch_group', glist)
        self.remove_zero_and_mask_taxon(self.option('welch_input'), './welch_group')
        two_group_test(self.work_dir + "/tmp_mask_otu.xls", "./welch_group",
        # two_group_test(self.option('welch_input').prop['path'], './welch_group',
                       self.work_dir + '/welch_result.xls', self.work_dir + '/welch_boxfile.xls', "welch",
                       str(1 - self.option('welch_ci')), self.option('welch_type'),
                       self.option('welch_correction'))
        self.logger.info(Config().SOFTWARE_DIR)
        cmd = self.r_path + " run_welch_test.r"
        self.logger.info("开始运行welch_T检验")
        command = self.add_command("welch_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("welch_cmd运行完成，开始运行计算置信区间")
            welch(self.work_dir + '/welch_result.xls', './welch_group', self.option('welch_coverage'))
            self.logger.info("生成单物种柱状图的数据")
            group_bar(self.work_dir + '/tmp_mask_otu.xls', './welch_group', self.work_dir + '/tmp_group_bar.xls', 'welch')
            # group_bar(self.option('welch_input').prop['path'], './welch_group', self.work_dir + '/welch_plot_group_bar.xls', 'welch')
            cmd1 = self.r_path + " run_welch_bar.r"
            bar_cmd = self.add_command("bar_cmd", cmd1).run()
            self.wait(bar_cmd)
            self.add_taxon(self.work_dir + '/tmp_group_bar.xls', self.work_dir + '/welch_plot_group_bar.xls')
            if bar_cmd.return_code == 0:
                self.logger.info("welch_cmd运行完成")
            else:
                self.set_error("bar_cmd运行出错!", code="32502907")
        else:
            self.set_error("welch_cmd运行出错!", code="32502908")

    def run_mann(self):
        # self.name_to_name = mask_taxon(self.option("mann_input").prop['path'], self.work_dir + "/tmp_mask_otu.xls")
        glist = [self.option('mann_gname')]
        self.option('mann_group').sub_group('./mann_group', glist)
        self.remove_zero_and_mask_taxon(self.option('mann_input'), './mann_group')
        two_group_test(self.work_dir + "/tmp_mask_otu.xls", "./mann_group",
        # two_group_test(self.option('mann_input').prop['path'], './mann_group',
                       self.work_dir + '/mann_result.xls', self.work_dir + '/mann_boxfile.xls', "mann",
                       str(1 - self.option('mann_ci')), self.option('mann_type'),
                       self.option('mann_correction'))
        cmd = self.r_path + " run_mann_test.r"
        self.logger.info("开始运行mann检验")
        command = self.add_command("mann_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            # self.logger.info("mann_cmd运行完成，开始运行计算置信区间")
            # student(self.work_dir + '/mann_result.xls', './mann_group', self.option('mann_coverage'))
            # bootstrap(self.option('mann_input').prop['path'], './mann_group', self.option('mann_coverage'))
            self.logger.info("生成单物种柱状图的数据")
            group_bar(self.work_dir + '/tmp_mask_otu.xls', './mann_group', self.work_dir + '/tmp_group_bar.xls', 'mann')
            # group_bar(self.option('mann_input').prop['path'], './mann_group', self.work_dir + '/mann_plot_group_bar.xls', 'mann')
            cmd1 = self.r_path + " run_mann_bar.r"
            bar_cmd = self.add_command("bar_cmd", cmd1).run()
            self.wait(bar_cmd)
            self.add_taxon(self.work_dir + '/tmp_group_bar.xls', self.work_dir + '/mann_plot_group_bar.xls')
            self.logger.info("开始运行计算置信区间")
            bootstrap(self.work_dir + '/mann_plot_group_bar.xls', './mann_group', self.option('mann_coverage'))
            if bar_cmd.return_code == 0:
                self.logger.info("mann_test运行完成")
            else:
                self.set_error("bar_cmd运行出错!", code="32502909")
        else:
            self.set_error("mann_cmd运行出错!", code="32502910")

    def run_signal(self):
        # self.name_to_name = mask_taxon(self.option("signal_input").prop['path'], self.work_dir + "/tmp_mask_otu.xls")
        glist = [self.option('signal_gname')]
        self.option('signal_group').sub_group('./signal_group', glist)
        self.remove_zero_and_mask_taxon(self.option('signal_input'), './signal_group')
        # two_group_test(self.option('signal_input').prop['path'], './signal_group',
        two_group_test(self.work_dir + '/tmp_mask_otu.xls', './signal_group',
                       self.work_dir + '/signal_result.xls', self.work_dir + '/signal_boxfile.xls', "signal",
                       str(1 - self.option('signal_ci')), self.option('signal_type'),
                       self.option('signal_correction'))
        cmd = self.r_path + " run_signal_test.r"
        self.logger.info("开始运行signal检验")
        command = self.add_command("signal_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            #self.logger.info("signal_cmd运行完成，开始运行计算置信区间")
            #bootstrap(self.option('signal_input').prop['path'], './signal_group', self.option('signal_coverage'))
            self.logger.info("生成单物种柱状图的数据")
            group_bar(self.work_dir + '/tmp_mask_otu.xls', './signal_group', self.work_dir + '/tmp_group_bar.xls', 'signal')
            # group_bar(self.option('signal_input').prop['path'], './signal_group', self.work_dir + '/signal_plot_group_bar.xls', 'signal')
            cmd1 = self.r_path + " run_signal_bar.r"
            bar_cmd = self.add_command("bar_cmd", cmd1).run()
            self.wait(bar_cmd)
            self.logger.info("开始运行计算置信区间")
            self.add_taxon(self.work_dir + '/tmp_group_bar.xls', self.work_dir + '/signal_plot_group_bar.xls')
            bootstrap(self.work_dir + '/signal_plot_group_bar.xls', './signal_group', self.option('signal_coverage'))
            if bar_cmd.return_code == 0:
                self.logger.info("signal_test运行完成")
            else:
                self.set_error("bar_cmd运行出错!", code="32502911")
        else:
            self.set_error("signal_cmd运行出错!", code="32502912")

    def run_kru(self):
        # self.name_to_name = mask_taxon(self.option("kru_H_input").prop['path'], self.work_dir + "/tmp_mask_otu.xls")
        glist = [self.option('kru_H_gname')]
        self.option('kru_H_group').sub_group('./kru_H_group', glist)
        self.remove_zero_and_mask_taxon(self.option('kru_H_input'), './kru_H_group')
        # mul_group_test(self.option('kru_H_input').prop['path'], self.work_dir + '/kru_H_result.xls',
        mul_group_test(self.work_dir + "/tmp_mask_otu.xls", self.work_dir + '/kru_H_result.xls',
                       self.work_dir + '/kru_H_boxfile.xls', './kru_H_group', "kru_H",
                       self.option('kru_H_correction'))
        cmd = self.r_path + " run_kru_H_test.r"
        self.logger.info("开始运行kru_H检验")
        command = self.add_command("kru_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("kru_cmd运行完成，开始运行post-hoc检验")
            self.posthoc(self.option("kru_H_methor"), self.work_dir + '/kru_H_result.xls', './kru_H_group',
                         self.option("kru_H_coverage"), './kru_H')
            self.logger.info("生成单物种柱状图的数据")
            group_bar(self.work_dir + '/tmp_mask_otu.xls', './kru_H_group', self.work_dir + '/tmp_group_bar.xls', 'kru_H')
            # group_bar(self.option('kru_H_input').prop['path'], './kru_H_group', self.work_dir + '/kru_H_plot_group_bar.xls', 'kru_H')
            cmd1 = self.r_path + " run_kru_H_bar.r"
            bar_cmd = self.add_command("bar_cmd", cmd1).run()
            self.wait(bar_cmd)
            self.add_taxon(self.work_dir + '/tmp_group_bar.xls', self.work_dir + '/kru_H_plot_group_bar.xls')
            if bar_cmd.return_code == 0:
                self.logger.info("kru_H_test运行完成")
            else:
                self.set_error("bar_cmd运行出错!", code="32502913")
        else:
            self.set_error("kru_cmd运行出错!", code="32502914")

    def run_anova(self):
        # self.name_to_name = mask_taxon(self.option("anova_input").prop['path'], self.work_dir + "/tmp_mask_otu.xls")
        glist = [self.option('anova_gname')]
        self.option('anova_group').sub_group('./anova_group', glist)
        self.remove_zero_and_mask_taxon(self.option('anova_input'), './anova_group')
        # mul_group_test(self.option('anova_input').prop['path'], self.work_dir + '/anova_result.xls',
        mul_group_test(self.work_dir + "/tmp_mask_otu.xls", self.work_dir + '/anova_result.xls',
                       self.work_dir + '/anova_boxfile.xls', './anova_group', "anova",
                       self.option('anova_correction'))
        cmd = self.r_path + " run_anova_test.r"
        self.logger.info("开始运行anova检验")
        command = self.add_command("anova_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("anova_cmd运行完成,开始运行post-hoc检验")
            self.posthoc(self.option("anova_methor"), self.work_dir + '/anova_result.xls',
                         './anova_group', self.option("anova_coverage"), './anova')
            self.logger.info("生成单物种柱状图的数据")
            group_bar(self.work_dir + '/tmp_mask_otu.xls', './anova_group', self.work_dir + '/tmp_group_bar.xls', 'anova')
            # group_bar(self.option('anova_input').prop['path'], './anova_group', self.work_dir + '/anova_plot_group_bar.xls', 'anova')
            cmd1 = self.r_path + " run_anova_bar.r"
            bar_cmd = self.add_command("bar_cmd", cmd1).run()
            self.wait(bar_cmd)
            self.add_taxon(self.work_dir + '/tmp_group_bar.xls', self.work_dir + '/anova_plot_group_bar.xls')
            if bar_cmd.return_code == 0:
                self.logger.info("anova_test运行完成")
            else:
                self.set_error("bar_cmd运行出错!", code="32502915")
        else:
            self.set_error("anova_cmd运行出错!", code="32502916")

    def posthoc(self, methor, statfile, groupfile, coverage, outfile):
        if methor == 'tukeykramer':
            tukeykramer(statfile, groupfile, coverage, outfile)
        if methor == 'gameshowell':
            gameshowell(statfile, groupfile, coverage, outfile)
        if methor == 'welchuncorrected':
            welchuncorrected(statfile, groupfile, coverage, outfile)
        if methor == 'scheffe':
            scheffe(statfile, groupfile, coverage, outfile)

    def replace_na(self, result_file):
        result_pd = pd.read_table(result_file, header=0)
        result_pd.fillna(1, inplace=True)
        output_file = os.path.join(self.output_dir, os.path.basename(result_file).replace('_new', ''))
        result_pd.to_csv(output_file, header=True, index=True, sep='\t')

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        for t in self.option('test').split(','):
            if t == 'chi':
                try:
                    self.add_taxon(self.work_dir + '/chi_result.xls', self.work_dir + '/chi_result_new.xls')
                    self.replace_na(self.work_dir + '/chi_result_new.xls')
                    # os.link(self.work_dir + '/chi_result_new.xls', self.output_dir + '/chi_result.xls')
                    self.add_taxon(self.work_dir + '/chi_CI.xls', self.work_dir + '/chi_CI_new.xls')
                    os.link(self.work_dir + '/chi_CI_new.xls', self.output_dir + '/chi_CI.xls')
                    # os.link(self.work_dir + '/chi_result.xls', self.output_dir + '/chi_result.xls')
                    # os.link(self.work_dir + '/chi_CI.xls', self.output_dir + '/chi_CI.xls')
                    self.logger.info("设置chi分析的结果目录成功")
                except:
                    self.logger.info("设置chi分析结果目录失败")
            elif t == 'fisher':
                try:
                    self.add_taxon(self.work_dir + '/fisher_result.xls', self.work_dir + '/fisher_result_new.xls')
                    self.replace_na(self.work_dir + '/fisher_result_new.xls')
                    # os.link(self.work_dir + '/fisher_result_new.xls', self.output_dir + '/fisher_result.xls')
                    self.add_taxon(self.work_dir + '/fisher_CI.xls', self.work_dir + '/fisher_CI_new.xls')
                    os.link(self.work_dir + '/fisher_CI_new.xls', self.output_dir + '/fisher_CI.xls')
                    # os.link(self.work_dir + '/fisher_result.xls', self.output_dir + '/fisher_result.xls')
                    # os.link(self.work_dir + '/fisher_CI.xls', self.output_dir + '/fisher_CI.xls')
                    self.logger.info("设置fisher分析的结果目录成功")
                except:
                    self.logger.info("设置fisher分析结果目录失败")
            elif t == 'student':
                try:
                    self.add_taxon(self.work_dir + '/student_result.xls', self.work_dir + '/student_result_new.xls')
                    self.replace_na(self.work_dir + '/student_result_new.xls')
                    # os.link(self.work_dir + '/student_result_new.xls', self.output_dir + '/student_result.xls')
                    self.add_taxon(self.work_dir + '/student_boxfile.xls', self.work_dir + '/student_boxfile_new.xls')
                    os.link(self.work_dir + '/student_boxfile_new.xls', self.output_dir + '/student_boxfile.xls')
                    self.add_taxon(self.work_dir + '/student_CI.xls', self.work_dir + '/student_CI_new.xls')
                    os.link(self.work_dir + '/student_CI_new.xls', self.output_dir + '/student_CI.xls')
                    # os.link(self.work_dir + '/student_result.xls', self.output_dir + '/student_result.xls')
                    # os.link(self.work_dir + '/student_boxfile.xls', self.output_dir + '/student_boxfile.xls')
                    # os.link(self.work_dir + '/student_CI.xls', self.output_dir + '/student_CI.xls')
                    self.logger.info("设置student分析的结果目录成功")
                except:
                    self.logger.info("设置student分析结果目录失败")
            elif t == 'welch':
                try:
                    self.add_taxon(self.work_dir + '/welch_result.xls', self.work_dir + '/welch_result_new.xls')
                    self.replace_na(self.work_dir + '/welch_result_new.xls')
                    # os.link(self.work_dir + '/welch_result_new.xls', self.output_dir + '/welch_result.xls')
                    self.add_taxon(self.work_dir + '/welch_boxfile.xls', self.work_dir + '/welch_boxfile_new.xls')
                    os.link(self.work_dir + '/welch_boxfile_new.xls', self.output_dir + '/welch_boxfile.xls')
                    self.add_taxon(self.work_dir + '/welch_CI.xls', self.work_dir + '/welch_CI_new.xls')
                    os.link(self.work_dir + '/welch_CI_new.xls', self.output_dir + '/welch_CI.xls')
                    # os.link(self.work_dir + '/welch_result.xls', self.output_dir + '/welch_result.xls')
                    # os.link(self.work_dir + '/welch_boxfile.xls', self.output_dir + '/welch_boxfile.xls')
                    # os.link(self.work_dir + '/welch_CI.xls', self.output_dir + '/welch_CI.xls')
                    self.logger.info("设置welch分析的结果目录成功")
                except:
                    self.logger.info("设置welch分析结果目录失败")
            elif t == 'mann':
                try:
                    self.add_taxon(self.work_dir + '/mann_result.xls', self.work_dir + '/mann_result_new.xls')
                    self.replace_na(self.work_dir + '/mann_result_new.xls')
                    # os.link(self.work_dir + '/mann_result_new.xls', self.output_dir + '/mann_result.xls')
                    self.add_taxon(self.work_dir + '/mann_boxfile.xls', self.work_dir + '/mann_boxfile_new.xls')
                    os.link(self.work_dir + '/mann_boxfile_new.xls', self.output_dir + '/mann_boxfile.xls')
                    # self.add_taxon(self.work_dir + '/bootstrap_CI.xls', self.work_dir + '/bootstrap_CI_new.xls')
                    # os.link(self.work_dir + '/bootstrap_CI_new.xls', self.output_dir + '/mann_CI.xls')
                    # os.link(self.work_dir + '/mann_result.xls', self.output_dir + '/mann_result.xls')
                    # os.link(self.work_dir + '/mann_boxfile.xls', self.output_dir + '/mann_boxfile.xls')
                    os.link(self.work_dir + '/bootstrap_CI.xls', self.output_dir + '/mann_CI.xls')
                    self.logger.info("设置mann分析的结果目录成功")
                except:
                    self.logger.info("设置mann分析结果目录失败")
            elif t == 'signal':
                try:
                    self.add_taxon(self.work_dir + '/signal_result.xls', self.work_dir + '/signal_result_new.xls')
                    self.replace_na(self.work_dir + '/signal_result_new.xls')
                    # os.link(self.work_dir + '/signal_result_new.xls', self.output_dir + '/signal_result.xls')
                    self.add_taxon(self.work_dir + '/signal_boxfile.xls', self.work_dir + '/signal_boxfile_new.xls')
                    os.link(self.work_dir + '/signal_boxfile_new.xls', self.output_dir + '/signal_boxfile.xls')
                    # self.add_taxon(self.work_dir + '/bootstrap_CI.xls', self.work_dir + '/bootstrap_CI_new.xls')
                    # os.link(self.work_dir + '/bootstrap_CI_new.xls', self.output_dir + '/signal_CI.xls')
                    # os.link(self.work_dir + '/signal_result.xls', self.output_dir + '/signal_result.xls')
                    # os.link(self.work_dir + '/signal_boxfile.xls', self.output_dir + '/signal_boxfile.xls')
                    os.link(self.work_dir + '/bootstrap_CI.xls', self.output_dir + '/signal_CI.xls')
                    self.logger.info("设置signal分析的结果目录成功")
                except:
                    self.logger.info("设置signal分析结果目录失败")
            elif t == 'anova':
                try:
                    self.add_taxon(self.work_dir + '/anova_result.xls', self.work_dir + '/anova_result_new.xls')
                    self.add_taxon(self.work_dir + '/anova_boxfile.xls', self.work_dir + '/anova_boxfile_new.xls')
                    self.replace_na(self.work_dir + '/anova_result_new.xls')
                    # os.link(self.work_dir + '/anova_result_new.xls', self.output_dir + '/anova_result.xls')
                    os.link(self.work_dir + '/anova_boxfile_new.xls', self.output_dir + '/anova_boxfile.xls')
                    # os.link(self.work_dir + '/anova_result.xls', self.output_dir + '/anova_result.xls')
                    # os.link(self.work_dir + '/anova_boxfile.xls', self.output_dir + '/anova_boxfile.xls')
                    for r, d, f in os.walk(self.work_dir, topdown=False):
                        filelist = f
                        self.logger.info(filelist)
                    for i in filelist:
                        if re.match(r'^anova_%s' % self.option("anova_methor"), i):
                            self.add_taxon(self.work_dir + '/' + i, self.work_dir + '/new_' + i)
                            os.link(self.work_dir + '/new_' + i, self.output_dir + '/' + i)
                            # os.link(self.work_dir + '/' + i, self.output_dir + '/' + i)
                    self.logger.info("设置anova分析的结果目录成功")
                except:
                    self.logger.info("设置anova分析结果目录失败")
            elif t == 'kru_H':
                try:
                    self.add_taxon(self.work_dir + '/kru_H_result.xls', self.work_dir + '/kru_H_result_new.xls')
                    self.add_taxon(self.work_dir + '/kru_H_boxfile.xls', self.work_dir + '/kru_H_boxfile_new.xls')
                    self.replace_na(self.work_dir + '/kru_H_result_new.xls')
                    # os.link(self.work_dir + '/kru_H_result_new.xls', self.output_dir + '/kru_H_result.xls')
                    os.link(self.work_dir + '/kru_H_boxfile_new.xls', self.output_dir + '/kru_H_boxfile.xls')
                    # os.link(self.work_dir + '/kru_H_result.xls', self.output_dir + '/kru_H_result.xls')
                    # os.link(self.work_dir + '/kru_H_boxfile.xls', self.output_dir + '/kru_H_boxfile.xls')
                    for r, d, f in os.walk(self.work_dir, topdown=False):
                        filelist = f
                    for i in filelist:
                        if re.match(r'^kru_H_%s' % self.option("kru_H_methor"), i):
                            self.add_taxon(self.work_dir + '/' + i, self.work_dir + '/new_' + i)
                            os.link(self.work_dir + '/new_' + i, self.output_dir + '/' + i)
                            # os.link(self.work_dir + '/' + i, self.output_dir + '/' + i)
                    self.logger.info("设置kru_H分析的结果目录成功")
                except:
                    self.logger.info("设置kru_H分析结果目录失败")
            elif t == 'estimator':
                try:
                    for r, d, f in os.walk(self.work_dir, topdown=False):
                        filelist = f
                    for i in filelist:
                        if re.match(r'^est_result', i):
                            os.link(self.work_dir + '/' + i, self.output_dir + '/' + i)
                    self.logger.info("设置est分析的结果目录成功")
                except:
                    self.logger.info("设置est分析结果目录失败")
