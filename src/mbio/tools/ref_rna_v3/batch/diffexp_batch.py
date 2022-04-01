# -*- coding: utf-8 -*-
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import unittest
from mbio.packages.ref_rna_v2.chart import Chart
from collections import defaultdict
import shutil
import json
__author__ = 'gdq'


class DiffexpBatchAgent(Agent):
    """
    Differential analysis based on edgeR, DEGseq, DEseq2.
    pvalue adjust way Integer in [1,2,3,4].
     1. Bonferroni. ---> Bonferroni;
     2. Bonferroni Step-down(Holm) ---> Holm;
     3. Benjamini and Hochberg False Discovery Rate ---> BH;
     4. FDR Benjamini-Yekutieli --->BY
     Default: 3
    """
    def __init__(self, parent):
        super(DiffexpBatchAgent, self).__init__(parent)
        options = [
            # count 和 exp 均为定量结果，其中count用于差异分析。
            dict(name="is_workflow", type="string", default="no"),
            dict(name="count", type="infile", format="ref_rna_v2.express_matrix"),
            dict(name="exp", type="infile", format="ref_rna_v2.express_matrix"),
            dict(name="exp_type", type="string", default="tpm"),

            # group为样本分组信息文件，要求至少两列，第一列为样本名，其他列为组名。
            # 没有重复实验时，可以用样本名作为组名。
            dict(name="group", type="infile", format="sample.group_table"),
            # cmp为比较信息文件，仅两列，第一列为对照组名，第二列为实验组名，要求与group信息保持一致。
            dict(name="cmp", type="infile", format="sample.control_table"),
            # dict(name="cmp", type="infile", format="denovo_rna_v2.compare_table"),
            # count_cutoff 和 over_cutoff 是一起用来判定某个基因是否可以进入后续差异分析的的依据。
            # 对于一个基因，N个样本中至少有over_cutoff个样本的表达超过count_cutoff的才能进入后续分析。
            dict(name="count_cutoff", type="int", default=4),
            dict(name="over_cutoff", type="int", default=None),
            #tpm_filter_threshold 20190522修改
            dict(name="tpm_filter_threshold", type="float", default=0),

            # 当用edgeR进行单样本和单样本比较时，需要输入经验值dispersion
            dict(name="dispersion", type="float", default=0.1),

            # pool为进程池的size. 用于差异分析的平行计算
            dict(name='pool', type='int', default=6),

            # pvalue 统计判断阈值
            dict(name="pvalue", type="float", default=0.05),
            dict(name="pvalue_padjust", type="string", default="padjust"),
            dict(name="fc", type="float", default=2),
            dict(name="padjust_way", type='string', default="BH"),

            # method: DESeq2, edgeR, DEGseq
            dict(name="method", type="string", default="DESeq2"),

            # output
            dict(name="output", type="string", default=None),
            dict(name="diff_result", type="string"),
            dict(name="diff_list", type="string", ),
            dict(name="diff_summary", type="string", ),
            # batch or pair
            dict(name="is_batch", type="bool", default=False),
            dict(name="has_batch", type="bool", default=True),
            dict(name="batch_matrix", type="infile", format="ref_rna_v2.common"),
            dict(name="is_duplicate", type="bool", default=True),
            dict(name='prob', type='float', default=0.8)

        ]
        self.add_option(options)
        self._memory_increase_step = 30

    def check_options(self):
        if self.option("method").lower() not in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma', 'noiseq']:
            raise OptionError("Method is incorrect", code = "33704301")
        if self.option('method').lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:
            pvalue = float(self.option("pvalue"))
            if not (0 <= pvalue <= 1):
                raise OptionError("pvalue cutoff must be in range [0-1]", code="33704302")
        fc = float(self.option('fc'))
        dispersion = float(self.option("dispersion"))
        if not (0 <= dispersion <= 1):
            raise OptionError("Dispersion argument of edgeR must be in range [0-1]", code = "33704303")
        if not fc >= 0:
            raise OptionError("fold change cutoff must be positive", code = "33704304")
        groups = self.option("group").prop['group_dict'].keys()
        cmp_groups = [x for y in self.option("cmp").prop['cmp_list'] for x in y]
        diff_groups = set(cmp_groups) - set(groups)
        if diff_groups:
            raise OptionError('Groups:%s in Compare is not in Group info', variables = (diff_groups), code = "33704305")

    def set_resource(self):
        cmp_num = len(self.option("cmp").prop['cmp_list'])
        if cmp_num <= self.option("pool"):
            self.option("pool", cmp_num)
        self._cpu = self.option("pool") + 1
        self._memory = "{}G".format(self.option("pool")*10)

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "差异分析结果目录"]
        #     ])
        # result_dir.add_regexp_rules([
        #     [r"*_vs_*.xls", "xls", "差异分析结果总表"],
        #     [r"*.DE.list", "xls", "差异基因列表"],
        #     [r"*summary.xls", "xls", "差异统计表"],
        #     ])
        super(DiffexpBatchAgent, self).end()


class DiffexpBatchTool(Tool):
    """
    Differential analysis based on edgeR, DEGseq, DEseq2.
    """
    def __init__(self, config):
        super(DiffexpBatchTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'program/Python/bin/python'
        # self.diff_toolbox = self.config.PACKAGE_DIR + "/ref_rna_v3/batch/diff_toolbox_batch.py"
        self.diff_toolbox = self.config.PACKAGE_DIR + "/ref_rna_v3/batch/diff_toolbox_batch_new.py"
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

        group_dict = dict()
        with open(self.option('group').path, 'r') as group:
            for line in group.readlines():
                if line.startswith('#') or not line.strip():
                    continue
                else:
                    sample, group = line.strip().split('\t')
                    group_dict.setdefault(group, list())
                    group_dict[group].append(sample)
        len_list = list()
        for g in group_dict:
            len_list.append(str(len(group_dict[g])))
            len_list.sort()
        if int(len_list[0]) == int(len_list[-1]) == 1:
            self.is_duplicate = False
        else:
            self.is_duplicate = True
        self.run_times = 0

    def diffexp(self):
        self.run_times += 1
        self.logger.info("本次分析为第{}次差异分析".format(str(self.run_times)))
        self.logger.info("本次分析使用软件{}".format(self.option("method")))
        self.logger.info("本次分析fc阈值为{}".format(str(self.option("fc"))))
        if not self.option("method").lower() == "noiseq":
            self.logger.info("本次分析使用{}作为显著性判断".format(self.option("pvalue_padjust")))
            self.logger.info("本次分析显著性阈值为{}".format(str(self.option("pvalue"))))
        cmd = '{} {} '.format(self.python_path, self.diff_toolbox)
        cmd += '-count {} '.format(self.option("count").prop['path'])
        cmd += '-exp {} '.format(self.option("exp").prop['path'])
        cmd += '--exp_type {} '.format(self.option("exp_type").lower())
        cmd += '-method {} '.format(self.option("method"))
        cmd += '-group {} '.format(self.option("group").prop['path'])
        cmd += '-cmp {} '.format(self.option("cmp").prop['path'])
        if self.option("output") is None:
            self.option("output", self.work_dir)
        else:
            if not os.path.exists(self.option("output")):
                os.mkdir(self.option("output"))
        cmd += '-output {} '.format(self.option("output"))
        cmd += '-pool {} '.format(self.option("pool"))
        if self.option('method').lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:
            cmd += '-pvalue {} '.format(self.option("pvalue"))
            cmd += '-sig_type {} '.format(self.option("pvalue_padjust"))
            if self.option("padjust_way").lower() == "bh":
                adjust_way = 3
            elif self.option("padjust_way").lower() == "bonferroni":
                adjust_way = 1
            elif self.option("padjust_way").lower() == "holm":
                adjust_way = 2
            elif self.option("padjust_way").lower() == "by":
                adjust_way = 4
            else:
                adjust_way = 3
            cmd += '-padjust_way {} '.format(adjust_way)

        cmd += '-fc {} '.format(self.option("fc"))
        cmd += '--dispersion {} '.format(self.option("dispersion"))
        cmd += '--count_cutoff {} '.format(self.option("count_cutoff"))
        cmd += '--tpm_filter_threshold {} '.format(self.option("tpm_filter_threshold"))
        if self.option("over_cutoff") is not None:
           cmd += '--passed_number_cutoff {} '.format(self.option("over_cutoff"))
        if self.option('is_batch') == True and self.option('has_batch') == True:
            cmd += '--batch {} '.format(self.option('batch_matrix').prop['path'])
        if self.option('method').lower() == 'noiseq':
            cmd += '--is_duplicate {} '.format(self.is_duplicate)
            cmd += '--prob {} '.format(self.option('prob'))
        #cmd += '--tpm_filter_threshold {}'.format(self.option("tpm_filter_threshold"))
        # cmd_name = 'diffexp'
        cmd_name = 'diffexp' + str(self.run_times)
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} 运行成功".format(cmd_name))
            else:
                self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704306")
        else:
            self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704307")
        self.set_output()


    def check_diff_results(self):
        diff_summary = glob.glob(self.option("output") + '/*summary*.xls')
        diff_summary_pd = pd.read_table(diff_summary[0], header=0, sep="\t", index_col=0)
        diff_summary_pd_remove_no = diff_summary_pd[(True - diff_summary_pd['sum'].isin([0]))]
        if diff_summary_pd_remove_no.shape[0] == 0:
            #min_diff_num = 0
            total_diff_num = 0
        else:
            total_diff_num = max([int(i) for i in list(diff_summary_pd_remove_no.iloc[0,])])
        self.logger.info("本次分析所有对比组共获得最差异数目为：{}".format(str(total_diff_num)))
        if self.is_duplicate:
            default_order = [("DESeq2","padjust",2.0,0.05),("DESeq2","pvalue",2.0,0.05),("edgeR","padjust",2.0,0.05),
                             ("edgeR","pvalue",2.0,0.05),("Limma","padjust",2.0,0.05),("Limma","pvalue",2.0,0.05),
                             ("DESeq2", "padjust", 1.5,0.05), ("DESeq2", "pvalue",1.5,0.05), ("edgeR", "padjust", 1.5,0.05),
                             ("edgeR", "pvalue", 1.5,0.05), ("Limma", "padjust", 1.5,0.05), ("Limma", "pvalue", 1.5,0.05),
                             ("NOIseq", 0.8, 2.0)]
        else:
            default_order = [("DEGseq", "padjust", 2.0,0.001), ("DEGseq", "pvalue", 2.0,0.001), ("edgeR", "padjust", 2.0,0.05),
                             ("edgeR", "pvalue", 2.0,0.05), ("Limma", "padjust", 2.0,0.05), ("Limma", "pvalue", 2.0,0.05),
                             ("DEGseq", "padjust", 1.5,0.05), ("DEGseq", "pvalue", 1.5,0.05), ("edgeR", "padjust", 1.5,0.05),
                             ("edgeR", "pvalue", 1.5,0.05), ("Limma", "padjust", 1.5,0.05), ("Limma", "pvalue", 1.5,0.05),
                             ("NOIseq", 0.8, 2.0)
                             ]
        if  self.option("method").lower() == "noiseq":
            # params_dict = ("method","prob","fc")
            params_pk = ("NOIseq",self.option("prob"),self.option("fc"))
        else:
            # params_dict = ( "method","pvalue_padjust","fc","pvalue")
            params_pk = (self.option("method"),self.option("pvalue_padjust"),self.option("fc"),self.option("pvalue"))
        if total_diff_num < 100:
            self.logger.info("本次分析所有对比组获得差异基因不足,重新调整参数")
            self.logger.info("本次分析调整前参数为：{}".format(str(params_pk)))
            if params_pk not in default_order:
                new_params_pk = default_order[0]
                params_dict = ("method", "pvalue_padjust", "fc", "pvalue")
                for key,value in zip(params_dict,new_params_pk):
                    self.option(key,value)

            else:
                index = default_order.index(params_pk)
                if index == len(default_order) -1:
                    self.set_error("尝试所有参数默认参数搭配方案,均未得到超过100个差异基因,,建议修改fc和pvalue值")
                else:
                    new_params_pk =  default_order[index+1]
                    self.logger.info("本次分析调整后参数为：{}".format(str(new_params_pk)))
                    if new_params_pk[0] != "NOIseq":
                        params_dict = ( "method","pvalue_padjust","fc","pvalue")
                    else:
                        params_dict = ("method", "prob", "fc")
                    for key, value in zip(params_dict, new_params_pk):
                        self.option(key, value)
            return False

        else:
            return True

    def modify_option_file(self):
        self.diff_options = defaultdict()
        if self.option('method').lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma','noiseq']:
            for option in ["method", "pvalue_padjust", "pvalue", "padjust_way", "prob", "fc"]:
                self.diff_options[option] = self.option(option)
        with open(os.path.join(self.work_dir,"option_file"),"w") as f:
            json.dump(self.diff_options, f, sort_keys=True, indent=4)

    def rerun_diffexp(self):
        all_deleted_inter_files = glob.glob(self.option("output") + '/*_vs_*.*') + glob.glob(self.option("output") + '/*summary*.xls')
        for file in all_deleted_inter_files:
            if os.path.isdir(file):
                shutil.rmtree(file)
            if  os.path.isfile(file):
                os.remove(file)
        self.diffexp()


    def set_output(self):
        if self.option("is_workflow") == "yes":
            result_usable = self.check_diff_results()
            if not result_usable:
                self.rerun_diffexp()
            else:
                self.modify_option_file()
                self.modify_result_file()
        else:
            self.modify_result_file()

    def modify_result_file(self):
        if self.option("method").lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:
            diff_files = glob.glob(self.option("output") + '/*_vs_*.*.xls')
            # diff_list = glob.glob(self.option("output") + '/*.DE.list')
            diff_summary = glob.glob(self.option("output") + '/*summary*.xls')
            diff_summary_pd = pd.read_table(diff_summary[0], header=0, sep="\t", index_col=0)
            diff_summary_pd_remove_no = diff_summary_pd[(True - diff_summary_pd['sum'].isin([0]))]
            diff_summary_pd_remove_no.to_csv(diff_summary[0], sep="\t")
            all_files = diff_files + diff_summary
            for each in all_files:
                fname = os.path.basename(each)
                if fname.endswith("normalize.xls"):
                    continue
                if fname.endswith("rawnormal.xls"):
                    continue
                if fname.endswith("normFactor.xls"):
                    continue
                if fname.endswith("sizeFactor.xls"):
                    continue
                link = os.path.join(self.output_dir, fname)
                if os.path.exists(link):
                    os.remove(link)
                os.link(each, link)
            diff_files_stat = glob.glob(os.path.join(self.output_dir, '*.{}.xls'.format(self.option("method").lower())))
            df1 = pd.DataFrame()
            # f_head = "seq_id\t \t \t"
            for diff_file_stat in diff_files_stat:
                if os.path.basename(diff_file_stat).startswith("total_diff_stat"):
                    continue
                elif os.path.basename(diff_file_stat).startswith("diff_summary_limma"):
                    continue
                ctrl, test = os.path.basename(diff_file_stat).split('.{}.xls'.format(self.option("method").lower()))[
                    0].split('_vs_')
                fname = os.path.basename(diff_file_stat).split(".")[0]
                # f_head = f_head + fname + "\t*\t*\t*\t*\t*\t*\t"
                df_t = pd.read_table(diff_file_stat, index_col="seq_id")
                df_t.rename(columns={"log2fc": "{}_log2fc({}/{})".format(fname, test, ctrl),
                                     "fc": "{}_fc({}/{})".format(fname, test, ctrl),
                                     "pvalue": "{}_pvalue".format(fname), "padjust": "{}_padjust".format(fname),
                                     "significant": "{}_significant".format(fname),
                                     "regulate": "{}_regulate".format(fname)},
                            inplace=True)
                df_core = pd.DataFrame(df_t, columns=["{}_fc({}/{})".format(fname, test, ctrl),
                                                      "{}_log2fc({}/{})".format(fname, test, ctrl),
                                                      "{}_pvalue".format(fname), "{}_padjust".format(fname),
                                                      "{}_significant".format(fname), "{}_regulate".format(fname)])
                df1 = pd.concat([df1, df_core], axis=1)
            df1.index.set_names("seq_id", inplace=True)
            df1.to_csv(self.output_dir + "/total_diff_stat.{}.xls".format(self.option("method").lower()), sep="\t")
        else:
            diff_files = glob.glob(self.option("output") + '/*_vs_*.*.xls')
            # diff_list = glob.glob(self.option("output") + '/*.DE.list')
            diff_summary = glob.glob(self.option("output") + '/*summary*.xls')
            diff_summary_pd = pd.read_table(diff_summary[0], header=0, sep="\t", index_col=0)
            diff_summary_pd_remove_no = diff_summary_pd[(True - diff_summary_pd['sum'].isin([0]))]
            diff_summary_pd_remove_no.to_csv(diff_summary[0], sep="\t")
            all_files = diff_files + diff_summary
            for each in all_files:
                fname = os.path.basename(each)
                if fname.endswith("normalize.xls"):
                    continue
                if fname.endswith("normFactor.xls"):
                    continue
                if fname.endswith("rawnormal.xls"):
                    continue
                if fname.endswith("sizeFactor.xls"):
                    continue
                link = os.path.join(self.output_dir, fname)
                if os.path.exists(link):
                    os.remove(link)
                os.link(each, link)
            diff_files_stat = glob.glob(os.path.join(self.output_dir, '*.{}.xls'.format(self.option("method").lower())))
            df1 = pd.DataFrame()
            # f_head = "seq_id\t \t \t"
            for diff_file_stat in diff_files_stat:
                ctrl, test = os.path.basename(diff_file_stat).split('.{}.xls'.format(self.option("method").lower()))[
                    0].split('_vs_')
                fname = os.path.basename(diff_file_stat).split(".")[0]
                # f_head = f_head + fname + "\t*\t*\t*\t*\t*\t*\t"
                df_t = pd.read_table(diff_file_stat, index_col="seq_id", sep='\t')
                df_t.rename(columns={"log2fc": "{}_log2fc({}/{})".format(fname, test, ctrl),
                                     "fc": "{}_fc({}/{})".format(fname, test, ctrl),
                                     "{}_mean".format(ctrl): "{}_{}_mean".format(fname, ctrl),
                                     "{}_mean".format(test): "{}_{}_mean".format(fname, test),
                                     "prob": "{}_prob".format(fname), "theta": "{}_theta".format(fname),
                                     "D": '{}_D'.format(fname),
                                     "significant": "{}_significant".format(fname),
                                     "regulate": "{}_regulate".format(fname)},
                            inplace=True)
                df_core = pd.DataFrame(df_t,
                                       columns=["{}_{}_mean".format(fname, ctrl), "{}_{}_mean".format(fname, test),
                                                "{}_fc({}/{})".format(fname, test, ctrl),
                                                "{}_log2fc({}/{})".format(fname, test, ctrl),
                                                "{}_theta".format(fname),
                                                '{}_D'.format(fname),
                                                "{}_prob".format(fname),
                                                "{}_significant".format(fname), "{}_regulate".format(fname)])
                df1 = pd.concat([df1, df_core], axis=1)
            df1.index.set_names("seq_id", inplace=True)
            df1.to_csv(self.output_dir + "/total_diff_stat.{}.xls".format(self.option("method").lower()), sep="\t")
        self.chart()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        diff_exp_summary = glob.glob(self.output_dir + '/*summary*.xls')[0]
        chart.chart_diffexp_stat(diff_exp_summary)

        for cmps in self.option("cmp").prop['cmp_list']:
            diff_exp = glob.glob(self.output_dir + '/*{}_vs_{}*.xls'.format(cmps[0], cmps[1]))[0]

            chart.chart_diffexp_scatter(diff_exp, "{}_vs_{}".format(cmps[0], cmps[1]), soft=self.option("method"), pvalue_padjust=self.option("pvalue_padjust"))
        chart.to_pdf()

        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            if os.path.exists(self.output_dir + "/" + os.path.basename(p)):
                os.remove(self.output_dir + "/" + os.path.basename(p))
            os.link(p, self.output_dir + "/" + os.path.basename(p))
        self.end()

    def run(self):
        super(DiffexpBatchTool, self).run()
        self.diffexp()
        # self.set_output()
        # self.chart()
        # self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    test_dir = '/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/tools/denovo_rna_v2/test_files'

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "Diffexp" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v3.batch.diffexp_batch",
            "instant": False,
            "options": dict(
                # count = '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/denova_rna_v2/unigene.count.matrix.xls',
                # exp = '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/denova_rna_v2/exp_matrix',
                # method = "DESeq2",
                # group = "/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/denova_rna_v2/group",
                # cmp = "/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/denova_rna_v2/cmp_little",
                count='/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/test20200710/data/unigene.count.matrix.xls',
                exp='/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/test20200710/data/exp_matrix',
                method="DESeq2",
                group="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/test20200710/data/group",
                cmp="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/test20200710/data/cmp",
                # count='/mnt/ilustre/users/sanger-dev/workspace/20200714/DiffexpBatch_tsg_37974_3873_5742/known_seqs_count.matrix',
                # exp='/mnt/ilustre/users/sanger-dev/workspace/20200714/DiffexpBatch_tsg_37974_3873_5742/exp_matrix',
                # method="DESeq2",
                # group="/mnt/ilustre/users/sanger-dev/workspace/20200714/DiffexpBatch_tsg_37974_3873_5742/group",
                # cmp="/mnt/ilustre/users/sanger-dev/workspace/20200714/DiffexpBatch_tsg_37974_3873_5742/cmp",
                pool=14,
                output=None,
                count_cutoff=5,
                dispersion=0.05,
                padjust_way="BH",
                pvalue=0.05,
                fc=2,
                # is_batch=True,
                # has_batch=True,
                # batch_matrix='/mnt/ilustre/users/sanger-dev/workspace/20200529/Refrna_tsg_37468/remote_input/pair_table/condition.txt'
            )
        }
        # data['options']['method'] = 'edgeR'
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()
        #
        # data['id'] += '1'
        # data['options']['method'] = 'DESeq2'
        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.run()
        # #
        # data['id'] += '2'
        # data['options']['method'] = 'DEGseq'
        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.run()


if __name__ == '__main__':
    unittest.main()


