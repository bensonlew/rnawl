# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import pandas as pd
from mbio.packages.itraq_and_tmt.metastat import *
import math
import glob
import subprocess
from collections import OrderedDict
import json
import numpy as np


class DiffAgent(Agent):
    def __init__(self, parent):
        super(DiffAgent, self).__init__(parent)
        options = [
            {"name":"group", "type":"infile","format":"itraq_and_tmt.group_table"}, # 分组文件
            {"name": "ratio_exp", "type": "infile", "format":"itraq_and_tmt.ratio_exp"},  # fisher,卡方，t检验的输入文件
            {"name":"cmp", "type":"infile","format":"itraq_and_tmt.compare_table"}, # 比较分组文件
            {"name": "fc_up", "type": "float", "default": 1.2},
            {"name": "fc_down", "type": "float", "default": 0.83},
            {'name': 'sig_type', 'type': 'string', 'default': 'pvalue'},
            {"name": "pvalue", "type": "float", "default": 0.05},
            {'name': 'padjust_way', 'type': 'int', 'default': 3},
            {"name": "correct_method", "type": "string", "default":"two.sided"}, #　卡方检验的页面没有单双尾检验
            {"name": "mul_test", "type": "string", "default": "none"},# 统一用p.adjust来矫正p值
            # param mul_test: 多重检验方法选择，默认为none，包括: ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"]
            {"name": "method_type", "type": "string", "default": "student"},
            # 费舍尔检验的选择单尾或双尾检验
            {"name": "group_dict", "type": "string", "default": "none"},
        ]
        self.add_option(options)
        self.step.add_steps("diff")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.diff.start()
        self.step.update()

    def step_end(self):
        self.step.diff.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数

        :return:
        """
        if not self.option("ratio_exp"):
            raise OptionError("必须设置输入文件:蛋白的scaled文件", code = "32502501")


    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["search_db.xls", " ", "搜库信息统计文件"],
        ])
        super(DiffAgent, self).end()


class DiffTool(Tool):
    def __init__(self, config):
        super(DiffTool, self).__init__(config)
        self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.3/bin/Rscript"
        self.cmds = []

    def split_group(self, group, cmp):
        df_group = pd.read_table(group, header=0,sep="\t")
        df_control = pd.read_table(cmp, header=0, sep="\t")
        control = df_control['#control'].tolist()
        other = df_control['other'].tolist()
        con_other = zip(control, other)
        i = 0
        for con in con_other:
            df_selected_control = df_group[(df_group.group == con[0])]
            df_selected_other = df_group[(df_group.group == con[1])]
            concat_info = [df_selected_control, df_selected_other]
            df_selected = pd.concat(concat_info)
            df_selected.to_csv(self.work_dir + "/cmp_revise" + str(i), sep='\t', index=False)
            i += 1
        # with open(cmp) as f_cmp:
        #     f_cmp.readline()
        #     i = 0
        #     for line in f_cmp:
        #         i += 1
        #         s_cmp = line.strip().split("\t")
        #         df_selected = df_group[(df_group.group == s_cmp[0]) | (df_group.group == s_cmp[1])]
        #         df_selected.to_csv(self.work_dir + "/cmp_revise" + str(i), sep = '\t', index=False)
        ## 注释掉的是以前有问题的代码，根本没有考虑顺序问题，果然没做过项目就是不行

    def diff(self):
        files = glob.glob(r"*cmp_revise*")
        i = 0
        for cmp in files:
            i += 1
            two_group_test(self.option("ratio_exp").prop['path'], cmp,self.work_dir + '/result'+ str(i) + '.xls',
             choose_test=self.option('method_type'),test_type=self.option("correct_method"))
            os.rename(self.work_dir + "/run_{}_test.r".format(self.option('method_type')), self.work_dir +
                      "/run_{}_test".format(self.option('method_type')) + str(i) + ".r")
            cmd = self.r_path + " run_{}_test".format(self.option('method_type')) + str(i) + ".r" #这里有一个空格，作为运行命令
            self.cmds.append(cmd) # 在使用command跑的时候，不仅cmd的名字需要变化，add.command(# xx,cmd)这个xx也需要变化
        return self.cmds

    def run_fisher(self):
        files = glob.glob(r"*cmp_revise*")
        i = 0
        for cmp in files:
            i += 1
            df = pd.read_table(cmp, header=0, sep="\t")
            fisher_sample1, fisher_sample2 = df["#sample"].tolist()
            two_sample_test(self.option("ratio_exp").prop['path'], self.work_dir + '/result'+ str(i) + '.xls', "fisher",
                fisher_sample1, fisher_sample2, test_type=self.option('correct_method'))
            os.rename(self.work_dir + "/run_fisher_test.r", self.work_dir + "/run_fisher_test" + str(i) + ".r")
            cmd = self.r_path + " run_fisher_test" + str(i) + ".r" #这里有一个空格，作为运行命令
            self.cmds.append(cmd) # 在使用command跑的时候，不仅cmd的名字需要变化，add.command(# xx,cmd)这个xx也需要变化
        return self.cmds

    def run_chi(self):
        files = glob.glob(r"*cmp_revise*")
        i = 0
        for cmp in files:
            i += 1
            df = pd.read_table(cmp, header=0, sep="\t")
            chi_sample1, chi_sample2 = df["#sample"].tolist()
            two_sample_test(self.option("ratio_exp").prop['path'], self.work_dir + '/result'+ str(i) + '.xls',  "chi",
                            chi_sample1, chi_sample2)
            os.rename(self.work_dir + "/run_chi_test.r", self.work_dir + "/run_chi_test" + str(i) + ".r")
            cmd = self.r_path + " run_chi_test" + str(i) + ".r" #这里有一个空格，作为运行命令
            self.cmds.append(cmd) # 在使用command跑的时候，不仅cmd的名字需要变化，add.command(# xx,cmd)这个xx也需要变化
        return self.cmds

    def run_cmds(self, func):
        """
        循环运行命令行，传入命令行列表
        """
        self.logger.info("开始运行命令！")
        count = 0
        if func == "diff":
            cmds = self.diff()
        elif func == "run_chi":
            cmds = self.run_chi()
        elif func == "run_fisher":
            cmds = self.run_fisher()
        for cmd in cmds:
            self.logger.info(cmd)
            self.logger.info("1234567890")
            try:
                subprocess.check_call(cmd, shell=True)
                self.logger.info(cmd + " 运行完成！")
                count += 1
            except subprocess.CalledProcessError:
                import traceback
                print traceback.format_exc()
                self.logger.info('CMD:{}'.format(cmd))
                break
        if count == len(cmds):
            self.logger.info("全部命令运行完成！！！")

    def diff_treat(self, fc_up, fc_down):
        """
        :param diff_r_result:
        diff函数调用r脚本得到的一个结果表格，只有表达均值，p值等相关信息，需要自己处理得到倍数，上下调等信息，这里作为输入文件
        fc_up, fc_down是页面端的输入参数，这里也是作为表格最终展示的一个筛选条件
        :return:
        """
        def my_test(row):
            return float(row[2]) / row[1]

        def my_test1(row):
            fc = float(row[2]) / row[1]
            return math.log(fc, 2)

        def my_test2(row):
            fc = float(row[2]) / row[1]
            if fc >= self.option("fc_up"):
             return "up"
            elif fc <= self.option("fc_down"):
             return "down"
            else:
             return "no change"

        def my_test3(row, p1=self.option("pvalue")):
            fc = float(row[2]) / row[1]
            if self.option('sig_type') == 'pvalue':
                p = row[3]
            if self.option('sig_type') == 'padjust':
                p = row[4]
            # if p <= p1 and (fc >= self.option("fc_up") or fc <= self.option("fc_down")):
            if p <= p1:
             return 'yes'
            if p > p1:
             return 'no'

        def my_test4(row):
            p = row[3]
            return math.log(p, 10)

        files = glob.glob(r"*result*")
        i = 0
        for file in files:
            i += 1
            df = pd.read_table(file, header=0, sep="\t")
            # df.drop(df.columns[[2, 4, 6]], axis=1, inplace=True)
            df.drop(df.columns[[2, 4, 6]], axis=1, inplace=True)
            columns = [x.split("-mean")[0] for x in list(df.columns)[0:3]]
            df.rename(columns={df.columns[0]: "accession_id", df.columns[1]: columns[1],
               df.columns[2]: columns[2]},inplace=True)
            df['padjust'] = self.multtest_correct(df['pvalue'], method=self.option('padjust_way'))
            df['fc'] = df.apply(my_test, axis=1)
            df['log2fc'] = df.apply(my_test1, axis=1)
            df['regulate'] = df.apply(my_test2, axis=1)
            df['significant'] = df.apply(my_test3, axis=1)
            str_1 = list(df.columns)[2] + "|" + list(df.columns)[1]
            df['log10pvalue'] = df.apply(my_test4, axis=1)
            df['log10padjust'] = df['padjust'].apply(lambda x: math.log(x, 10))
            df['compare'] = str_1
            # df = df[(df.fc >= self.option("fc_up")) | (df.fc ==self.option("fc_down"))]
            df.to_csv(self.work_dir + "/diffcmp" + str(i) + ".csv", sep = '\t', index=False)

    def diff_treat_fisher_chi(self, fc_up, fc_down):
        """
        :param diff_r_result:
        diff函数调用r脚本得到的一个结果表格，只有表达均值，p值等相关信息，需要自己处理得到倍数，上下调等信息，这里作为输入文件
        fc_up, fc_down是页面端的输入参数，这里也是作为表格最终展示的一个筛选条件
        :return:
        """
        def my_test(row):
            return float(row[2]) / row[1]

        def my_test1(row):
            fc = float(row[2]) / row[1]
            return math.log(fc, 2)

        def my_test2(row):
            fc = float(row[2]) / row[1]
            if fc >= self.option("fc_up"):
             return "up"
            elif fc <= self.option("fc_down"):
             return "down"
            else:
             return "no change"

        def my_test3(row, p1=self.option("pvalue")):
            fc = float(row[2]) / row[1]
            if self.option('sig_type') == 'pvalue':
                p = row[3]
            if self.option('sig_type') == 'padjust':
                p = row[4]
            if p <= p1 :
             return 'yes'
            if p > p1:
             return 'no'

        def my_test4(row):
            p = row[3]
            return math.log(p, 10)

        try:
            group_dict = json.loads(self.option('group_dict'), object_pairs_hook=OrderedDict)
        except:
            group_dict = self.option('group').prop["group_dict"]
        new_dict = {str(v[0]):k for k,v in group_dict.items()}
        print(new_dict)
        print(999999999)
        files = glob.glob(r"*result*")
        i = 0
        for file in files:
            i += 1
            df = pd.read_table(file, header=0, sep="\t")
            columns = [new_dict[x.split("-propotion")[0]] for x in list(df.columns)[1:3]]
            df.rename(columns={df.columns[0]: "accession_id", df.columns[1]:columns[0],
               df.columns[2]: columns[1]},inplace=True)
            df.drop(df.columns[[4]], axis=1, inplace=True)
            df['padjust'] = self.multtest_correct(df['pvalue'], method=self.option('padjust_way'))
            df['fc'] = df.apply(my_test, axis=1)
            df['log2fc'] = df.apply(my_test1, axis=1)
            df['regulate'] = df.apply(my_test2, axis=1)
            df['significant'] = df.apply(my_test3, axis=1)
            str_1 = list(df.columns)[2] + "|" + list(df.columns)[1]
            df['log10pvalue'] = df.apply(my_test4, axis=1)
            df['log10padjust'] = df['padjust'].apply(lambda x: math.log(x, 10))
            df['compare'] = str_1
            # df = df[(df.fc >= self.option("fc_up")) | (df.fc ==self.option("fc_down"))]
            df.to_csv(self.work_dir + "/diffcmp" + str(i) + ".csv", sep = '\t', index=False)


    def merge(self):
        """
        根据accession_id合并所有的文字名里面带有diffcmp表格
        :return:
        也可以先读出来一个，在循环，不过需要用位置
        res = pd.read_table(files[0], header=0, sep='\t')
        for i in range(1,len(files)):
            res.merge(pd.read_table(f[i], header=0, sep='\t'), how='outer', on='accession_id')

        """
        files = glob.glob(r"*diffsummary*")
        i = 0
        list_df = list()
        for file in files:
            i += 1
            file_name = pd.read_table(file, header=0, sep="\t",
                                      index_col=0)# concat按照索引合并，即使你是乱的，也可以拼起来
            file_name.columns = [x.replace("_vs_", "|") for x in list(file_name.columns)]
            list_df.append(file_name)
        df_merged = pd.concat(list_df, axis=1)
        # df_merged = reduce(lambda  left,right: pd.merge(left,right, on=['accession_id'], how='outer'), list_df)
        df_merged.index.name = "accession_id"
        df_merged.to_csv(self.work_dir + "/allsummary" + ".xls", sep = '\t')

    def statistic_allsummary(self):
        with open(self.work_dir + "/allsummary.xls", "r") as s4:
            head_indel = s4.readline()
            head_indel_list = head_indel.strip("\n").split("\t")[1:]
            targets = [x.strip() for x in "yes, no".split(',')]
            data_dict = OrderedDict(zip(head_indel_list, [{x:0 for x in targets} for _ in head_indel_list]))

            for line in s4:
                tmp_list = line.strip("\n").split("\t")[1:] #一定要strip,
                # 要不然最后一个统计会出现全是0的情况
                for ind, sample in enumerate(tmp_list):
                    if sample == 'yes':
                        data_dict[head_indel_list[ind]]['yes'] += 1
                    elif sample == 'no':
                        data_dict[head_indel_list[ind]]['no'] += 1

            result = pd.DataFrame(data_dict)
            result.index.name = 'type'
            result.to_csv(self.work_dir + "/diff_statistic", sep = "\t")

    def get_cmp_dict(self):
        files = glob.glob(r"*cmp_revise*")
        i = 0
        list_df = list()
        dict_cmp_1 = OrderedDict()
        dict_final = OrderedDict()
        for file in files:
            i += 1
            df = pd.read_table(file, header=0, sep="\t")#
            # concat按照索引合并，即使你是乱的，也可以拼起来
            a = df.iloc[:, 1].tolist()
            group_1, group_2 = [x for i, x in enumerate(a) if a.index(x) == i]
            # group_1 = list(df.iloc[:, 1])[len(list(df.iloc[:, 0]))/2 - 1]
            # group_2 = list(df.iloc[:, 1])[len(list(df.iloc[:, 0]))/2 + 1]
            str1 = group_2 + "|" + group_1
            samples = list(df.iloc[:, 0])
            dict_cmp_1[str1] = samples
            dict_final.update(dict_cmp_1)
        # result = pd.DataFrame(dict_final)
        # 这样就可以避免字典的键的值长度不一样的情况
        result = pd.DataFrame(dict([(k, pd.Series(v)) for k,v in dict_final.items()]))
        result.index.name = 'cmp'
        result.to_csv(self.work_dir + "/compare_dict.xls", sep = "\t", index=False)

    def nosig_up_down(self):
        df = pd.read_table(self.work_dir + "/diff_statistic", sep = "\t", header=0)
        groups = list(df.columns)[1::3]
        ups = list(df.columns)[2::3]
        downs = list(df.columns)[3::3]
        nosig = [df.iloc[0, x] for x in range(1, df.shape[1], 3)]
        up_num = [df.iloc[1, x] for x in range(2, df.shape[1], 3)]
        down_num = [df.iloc[1, x] for x in range(3, df.shape[1], 3)]
        nosig_str = ["nosig_" + str(i) for i in nosig]
        up_str = ["up_" + str(i) for i in up_num]
        down_str = ["down_" + str(i) for i in down_num]
        zip_tuple = zip(groups, nosig_str, down_str, up_str)
        dict_zip = OrderedDict()
        for i in zip_tuple:
            dict_zip[i[0]] = i[1:]
        result = pd.DataFrame(dict_zip)
        result.to_csv(self.work_dir + "/num_summary.xls", sep = "\t", index=False)

    def diff_summary(self):
        def my_test(row):
            if row[8] == "yes" and (row[7] == "up" or row[7] == "down"):
                return "yes"
            else:
                return "no"

        def my_test1(row):
            if row[7] == "up" and row[8] == "yes":
                return "yes"
            if row[7] == "up" and row[8] == "no":
                return "no"

        def my_test2(row):
            if row[7] == "down" and row[8] == "yes":
                return "yes"
            if row[7] == "down" and row[8] == "no":
                return "no"

        files = glob.glob(r"*diffcmp*")
        i = 0
        for file in files:
            i += 1
            df = pd.read_table(file, header=0, sep="\t")
            columns = list(df.columns)[1:3]
            self.logger.info(columns)
            self.logger.info("双击8888888")
            str1 = columns[1] + "_" + "vs" + "_" + columns[0]
            str1_up = str1 + "_up"
            str1_down = str1 + "_down"
            df[str1] = df.apply(my_test, axis=1)
            df[str1_up] = df.apply(my_test1, axis=1)
            df[str1_down] = df.apply(my_test2, axis=1)
            df_new = df.loc[:, ["accession_id", str1, str1_up, str1_down]]
            df_new.fillna("no", inplace=True)
            df_new.index.name = "accession_id"
            # df = df[(df.fc >= self.option("fc_up")) | (df.fc ==self.option("fc_down"))]
            df_new.to_csv(self.work_dir + "/diffsummary" + str(i) + ".csv",
                          sep = '\t', index=False)

    def multtest_correct(self, p_values, method=3):
        """
        1. Bonferroni. ---> bonferroni
        2. Bonferroni Step-down(Holm) ---> Holm
        3. Benjamini and Hochberg False Discovery Rate ---> BH
        4. FDR Benjamini-Yekutieli --->BY
        :param pvalue_list:
        :param method:
        :return: np.array
        """
        pvalue_list = list(p_values)
        n = len(pvalue_list)
        if method == 1:
            fdr = [eachP*n for eachP in pvalue_list]
        elif method == 2:
            sorted_pvalues = sorted(pvalue_list)
            fdr = [eachP*(n - sorted_pvalues.index(eachP)) for eachP in pvalue_list]
        elif method == 3:
            sorted_pvalues = sorted(pvalue_list)
            fdr = [eachP*n/(sorted_pvalues.index(eachP)+1) for eachP in pvalue_list]
        elif method == 4:
            _, fdr = self.fdr_correction(pvalue_list, alpha=0.05, method='negcorr', is_sorted=False)
        fdr = np.array(fdr)
        fdr[fdr > 1] = 1.
        return fdr

    def fdr_correction(self, pvals, alpha=0.05, method='indep', is_sorted=False):
        """pvalue correction for false discovery rate
        This covers Benjamini/Hochberg for independent or positively correlated and
        Benjamini/Yekutieli for general or negatively correlated tests. Both are
        available in the function multipletests, as method=`fdr_bh`, resp. `fdr_by`.
        Parameters
        ----------
        pvals : array_like
            set of p-values of the individual tests.
        alpha : float
            error rate
        method : {'indep', 'negcorr')
        Returns
        -------
        rejected : array, bool
            True if a hypothesis is rejected, False if not
        pvalue-corrected : array
            pvalues adjusted for multiple hypothesis testing to limit FDR
        Notes
        -----
        If there is prior information on the fraction of true hypothesis, then alpha
        should be set to alpha * m/m_0 where m is the number of tests,
        given by the p-values, and m_0 is an estimate of the true hypothesis.
        (see Benjamini, Krieger and Yekuteli)
        """
        def _ecdf(x):
            """
            no frills empirical cdf used in fdrcorrection
            """
            nobs = len(x)
            return np.arange(1, nobs+1)/float(nobs)

        pvals = np.asarray(pvals)
        if not is_sorted:
            pvals_sortind = np.argsort(pvals)
            pvals_sorted = np.take(pvals, pvals_sortind)
        else:
            pvals_sorted = pvals  # alias

        if method in ['i', 'indep', 'p', 'poscorr']:
            ecdffactor = _ecdf(pvals_sorted)
        elif method in ['n', 'negcorr']:
            cm = np.sum(1./np.arange(1, len(pvals_sorted)+1))   #corrected this
            ecdffactor = _ecdf(pvals_sorted) / cm
        else:
            raise ValueError('only indep and negcorr implemented')
        reject = pvals_sorted <= ecdffactor*alpha
        if reject.any():
            rejectmax = max(np.nonzero(reject)[0])
            reject[:rejectmax] = True

        pvals_corrected_raw = pvals_sorted / ecdffactor
        pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
        del pvals_corrected_raw
        pvals_corrected[pvals_corrected > 1] = 1
        if not is_sorted:
            pvals_corrected_ = np.empty_like(pvals_corrected)
            pvals_corrected_[pvals_sortind] = pvals_corrected
            del pvals_corrected
            reject_ = np.empty_like(reject)
            reject_[pvals_sortind] = reject
            return reject_, pvals_corrected_
        else:
            return reject, pvals_corrected

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        diff_summary = self.work_dir + "/" + "allsummary.xls"
        if os.path.exists(os.path.join(self.output_dir, "diff_protein_summary.xls")):
            os.remove(os.path.join(self.output_dir, "diff_protein_summary.xls"))
        os.link(diff_summary, os.path.join(self.output_dir, "diff_protein_summary.xls"))
        for f in os.listdir(self.work_dir):
            if f.startswith("diffcmp"):
                with open(f, 'r') as f1:
                    f1.readline()
                    lines = f1.readline()
                    line = lines.strip().split("\t")
                    cmp = line[-1].split("|")
                    compare = '{}_vs_{}'.format(cmp[0], cmp[1])
                    if os.path.exists(os.path.join(self.output_dir, compare + "_diff.xls")):
                        os.remove(os.path.join(self.output_dir, compare + "_diff.xls"))
                    os.link(self.work_dir + "/" + f, os.path.join(self.output_dir, compare + "_diff.xls"))
        self.logger.info("设置搜库结果目录")

    def run(self):
        super(DiffTool, self).run()
        self.split_group(self.option("group").prop['path'], self.option("cmp").prop['path'])
        if self.option('method_type') == 'student' or self.option(
                'method_type') == 'welch' or self.option('method_type') == 'signal':
            # signal符号检验用于配对t检验
            func = "diff"
            self.run_cmds(func)
            self.diff_treat(self.option("fc_up"), self.option("fc_down"))
        elif self.option('method_type') == 'chi':
            func = "run_chi"
            self.run_cmds(func)
            self.diff_treat_fisher_chi(self.option("fc_up"), self.option("fc_down"))
        elif self.option('method_type') == 'fisher':
            func = "run_fisher"
            self.run_cmds(func)
            self.diff_treat_fisher_chi(self.option("fc_up"), self.option("fc_down"))

        self.diff_summary()
        self.merge()
        self.statistic_allsummary()
        self.get_cmp_dict()
        self.nosig_up_down()
        self.set_output()
        self.end()
