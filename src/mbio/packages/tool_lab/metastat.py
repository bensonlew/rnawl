#!/mnt/ilustre/users/sanger/app/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "qiuping"


from mako.template import Template
from biocluster.config import Config
import os


this_file_dir = os.path.dirname(os.path.realpath(__file__))
def two_group_test(inputfile, groupfile, outputfile, boxfile, choose_test, ci=0.95, test_type="two.side", mul_test="none", norm="T"):
    """
    生成并运行R脚本，进行两组样品的差异性分析，包括student T检验，welch T检验，wilcox秩和检验

    :param inputfile: 输入的某一水平的otu_taxon_table
    :param groupfile: 输入分组文件
    :param outputfile: 输出的结果文件
    :param choose_test：选择两组检验的分析方法，包括：["student","welch","mann"]
    :param ci: 置信区间水平，默认为0.95
    :param test_type: 选择单双尾检验，默认为two.side，包括：["two.side","less","greater"]
    :param mul_test: 多重检验方法选择，默认为none，包括: ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"]
    :param norm: 代谢新增参数，计算时不归一化, # add by shaohua.yuan 20180801
    """
    f = Template(filename= this_file_dir + '/two_group_test.r')
    two_test = f.render(inputfile=inputfile, outputfile=outputfile, boxfile=boxfile, groupfile=groupfile, choose_test=choose_test, test_type=test_type, mul_test=mul_test, ci=ci,norm=norm)
    with open("run_%s_test.r" % choose_test, 'w') as rfile:
        rfile.write("%s" % two_test)


def two_sample_test(inputfile,outputfile,choose_test, sample1, sample2, ci=0.95, test_type="two.side", mul_test="none", cpu=1):
    """

    生成并运行R脚本，进行两样品的卡方检验或者费舍尔检验
    :param inputfile: 输入的某一水平的otu_taxon_table
    :param sample1: 输入样本1
    :param sample2: 输入样本2
    :param outputfile: 输出的结果文件
    :param choose_test：选择两组检验的分析方法，包括：["chi_sq","fisher"]
    :param ci: 置信区间水平，默认为0.95
    :param test_type: 选择单双尾检验，默认为two.side，包括：["two.side","less","greater"]
    :param mul_test: 多重检验方法选择，默认为none，包括: ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"]
    :param cpu: 子进程 cpu数
    """
    f = Template(filename= this_file_dir + '/two_sample_test.r')
    two_test = f.render(inputfile=inputfile, outputfile=outputfile, choose_test=choose_test, sample1=sample1, sample2=sample2, ci=ci, test_type=test_type, mul_test=mul_test, cpu_num=cpu)
    with open("run_%s_test.r" % choose_test, 'w') as rfile:
        rfile.write('%s' % two_test)


def mul_group_test(inputfile, outputfile, boxfile, groupfile, choose_test, mul_test="none",need_percents="T",percent_abund="T"):
    """
    生成并运行R脚本，进行多组样本的差异性分析，包括克鲁斯卡尔-Wallis秩和检验、anova分析
    :param inputfile: 输入的某一水平的otu_taxon_table
    :param groupfile: 输入分组文件
    :param outputfile: 输出的结果文件
    :param choose_test：选择两组检验的分析方法，包括：["kru_H","anova"]
    :param mul_test: 多重检验方法选择，默认为none，包括: ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"]
    """
    f = Template(filename=this_file_dir + '/mul_group_test.r')
    mul_test = f.render(inputfile=inputfile, outputfile=outputfile, boxfile=boxfile, groupfile=groupfile, choose_test=choose_test, mul_test=mul_test,need_percents=need_percents,percent_abund=percent_abund)
    with open("run_%s_test.r" % choose_test, 'w') as rfile:
        rfile.write("%s" % mul_test)

def mul_group_test_toollab(inputfile, outputfile, boxfile, groupfile, choose_test, mul_test="none", need_percents="true", percent_abund="true"):
    """
    生成并运行R脚本，进行多组样本的差异性分析，包括克鲁斯卡尔-Wallis秩和检验、anova分析
    :param inputfile: 输入的某一水平的otu_taxon_table
    :param groupfile: 输入分组文件
    :param outputfile: 输出的结果文件
    :param choose_test：选择两组检验的分析方法，包括：["kru_H","anova"]
    :param mul_test: 多重检验方法选择，默认为none，包括: ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"]
    :param choose_test：结果是否需要百分比形式
    """
    f = Template(filename=this_file_dir + '/mul_group_test_toollab.r')
    mul_test = f.render(inputfile=inputfile, outputfile=outputfile, boxfile=boxfile, groupfile=groupfile, choose_test=choose_test, mul_test=mul_test, need_percents=need_percents, percent_abund=percent_abund)
    with open("run_%s_test.r" % choose_test, 'w') as rfile:
        rfile.write("%s" % mul_test)

def est_ttest(inputfile, outputfile, groupfile, choose_test):
    """
    生成并运行R脚本，进行alpha指数T检验分析
    :param inputfile: 输入的某一水平的多样性指数表
    :param groupfile: 输入分组文件
    :param outputfile: 输出的结果文件
    """
    f = Template(filename=this_file_dir + '/alpha_ttest.r')
    mul_test = f.render(inputfile=inputfile, outputfile=outputfile, groupfile=groupfile, choose_test=choose_test)
    with open("run_est_ttest.r", 'w') as rfile:
        rfile.write("%s" % mul_test)


def group_bar(inputfile, groupfile, outputfile, choose_test):
    '''
    计算单物种柱形图的数据
    '''
    f = Template(filename=this_file_dir + '/group_bar.r')
    mul_test = f.render(inputfile=inputfile, outputfile=outputfile, groupfile=groupfile)
    with open("run_{}_bar.r".format(choose_test), 'w') as rfile:
        rfile.write("%s" % mul_test)

def group_bar_toollab(inputfile, groupfile, outputfile, choose_test, method="T"):
    '''
    计算单物种柱形图的数据
    '''
    f = Template(filename=this_file_dir + '/group_bar_toollab.r')
    mul_test = f.render(inputfile=inputfile, outputfile=outputfile, groupfile=groupfile,method=method)
    with open("run_{}_bar.r".format(choose_test), 'w') as rfile:
        rfile.write("%s" % mul_test)

def diff_kegg_test(proteinfile, enzymefile, filterfile, groupfile, outputfile, choose_test, ci=0.95, test_type="two.side", mul_test="none"):
    """
    计算多个分组中两两组之间的差异检验
    :param proteinfile: 输入kegg蛋白水平的丰度
    :param enzymefile: 输入kegg酶水平的丰度
    :param filterfile: 用于筛选蛋白与酶的表格，第一行为蛋白，第二行为酶
    :param groupfile: 输入分组文件
    :param outputfile: 输出的结果文件
    :param choose_test: 选择两组检验的分析方法，包括：["student", "welch", "mann"]
    :param ci: 置信区间水平，默认为0.95
    :param test_type: 选择单双尾检验，默认为two.side，包括：["two.side", "less", "greater"]
    :param mul_test: 多重检验方法选择，默认为none，包括：["holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"]
    :param norm: 计算时是否进行归一化
    """
    f = Template(filename=this_file_dir + '/diff_kegg_test.r')
    kegg_test = f.render(proteinfile=proteinfile, enzymefile=enzymefile, filterfile=filterfile, groupfile=groupfile, outputfile=outputfile, choose_test=choose_test, ci=ci, test_type=test_type, mul_test=mul_test)
    with open("run_%s_test.r" % choose_test, 'w') as rfile:
        rfile.write("%s" % kegg_test)