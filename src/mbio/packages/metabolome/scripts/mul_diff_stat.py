# -*- coding: utf-8 -*-
# __author__ = "shaohua.yuan"
# last_modifiy = modified 2018.06.13

from mako.template import Template
from biocluster.config import Config
import os, sys, argparse
import pandas as pd


class MulDiffStat(object):
    """
    多元统计pca,plsda,oplsda分析
    """
    def __init__(self):
        self.this_file_dir = os.path.dirname(os.path.realpath(__file__))

    def mul_pca(self, inputfile, output, groupfile, ci="0.95", data_trans="standard", mul_type="pca", perm=""):
        """
        生成并运行R脚本，进行两组样品的差异性分析，包括student T检验，welch T检验，wilcox秩和检验

        :param inputfile: 输入的某一水平的otu_taxon_table
        :param groupfile: 输入分组文件
        :param outputfile: 输出的结果文件
        :param choose_test：选择两组检验的分析方法，包括：["student","welch","mann"]
        :param ci: 置信区间水平，默认为0.95
        :param test_type: 选择单双尾检验，默认为two.side，包括：["two.side","less","greater"]
        :param mul_test: 多重检验方法选择，默认为none，包括: ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"]
        """
        f = Template(filename=self.this_file_dir + '/mul_diff_stat.r')
        two_test = f.render(inputfile=inputfile, output=output, groupfile=groupfile, data_trans=data_trans,
                            mul_type=mul_type, ci=ci, perm=perm, scr_dir=self.this_file_dir)
        mul_type = mul_type.replace(";", "_")
        group = groupfile.split("/")[-1].rpartition("_group")[0]
        rname = "run_{}_{}.r".format(mul_type, group)
        with open(rname, 'w') as rfile:
            rfile.write("%s" % two_test)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, metavar="exp_matrix_file", required=True, help="expression matrix file")
    #parser.add_argument('-group', type=str, metavar="group_file", default=None, help="file with two col: sample\tgroup")
    parser.add_argument('-o', type=str, required=True, help="default is local dir. Output directory name")
    parser.add_argument('-g', type=str, required=True, metavar="group_file", help="group_file")
    parser.add_argument('-dg', metavar='diff_group_name', help="diff group name")
    parser.add_argument('-ci', metavar='confidence', type=str, default="0.95", help="confidence")
    parser.add_argument('-trans', metavar='scale_method', default='standard', help="data scale method")
    parser.add_argument('-mt', metavar='mul_type', default='pca', help="pca,plsda,oplsda")
    parser.add_argument('-perm', metavar='permutation', type=str, default="200", help="permutation number")
    args = parser.parse_args()
    exp_pd = args.i
    outDir = args.o
    group = args.g
    ci = args.ci
    data_trans = args.trans
    mul_type = args.mt
    perm = args.perm
    if args.dg:
        diff = args.dg
    run = MulDiffStat()
    run_pca = run.mul_pca(exp_pd, outDir, group, ci=ci, data_trans=data_trans, mul_type=mul_type, perm=perm)

