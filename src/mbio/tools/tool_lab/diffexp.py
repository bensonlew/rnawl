# -*- coding: utf-8 -*-
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import unittest
import re
import numpy as np
from mbio.packages.whole_transcriptome.utils import runcmd
__author__ = 'gdq'


class DiffexpAgent(Agent):
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
        super(DiffexpAgent, self).__init__(parent)
        options = [
            # count 和 exp 均为定量结果，其中count用于差异分析。
            dict(name="count", type="infile", format="ref_rna_v2.common"),
            # dict(name="exp", type="infile", format="ref_rna_v2.express_matrix"),
            # dict(name="exp_type", type="string", default="tpm"),

            # group为样本分组信息文件，要求至少两列，第一列为样本名，其他列为组名。
            # 没有重复实验时，可以用样本名作为组名。
            dict(name="group", type="infile", format="ref_rna_v2.common"),
            # cmp为比较信息文件，仅两列，第一列为对照组名，第二列为实验组名，要求与group信息保持一致。
            dict(name="cmp", type="infile", format="ref_rna_v2.common"),

            # # count_cutoff 和 over_cutoff 是一起用来判定某个基因是否可以进入后续差异分析的的依据。
            # # 对于一个基因，N个样本中至少有over_cutoff个样本的表达超过count_cutoff的才能进入后续分析。
            # dict(name="count_cutoff", type="int", default=4),
            # dict(name="over_cutoff", type="int", default=None),


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
            dict(name='prob', type='float', default=0.8),
            dict(name='norm', type='string'),
            dict(name='model', type='string'),
            dict(name='test', type='string'),
            dict(name='voom', type='bool'),

        ]
        self.add_option(options)

    def check_options(self):
        # if self.option("method").lower() not in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma', 'noiseq']:
        #     raise OptionError("Method is incorrect", code = "33704301")
        # if self.option('method').lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:
        #     pvalue = float(self.option("pvalue"))
        #     if not (0 <= pvalue <= 1):
        #         raise OptionError("pvalue cutoff must be in range [0-1]", code="33704302")
        # fc = float(self.option('fc'))
        # dispersion = float(self.option("dispersion"))
        # if not (0 <= dispersion <= 1):
        #     raise OptionError("Dispersion argument of edgeR must be in range [0-1]", code = "33704303")
        # if not fc >= 0:
        #     raise OptionError("fold change cutoff must be positive", code = "33704304")
        # groups = self.option("group").prop['group_dict'].keys()
        # cmp_groups = [x for y in self.option("cmp").prop['cmp_list'] for x in y]
        # diff_groups = set(cmp_groups) - set(groups)
        # if diff_groups:
        #     raise OptionError('Groups:%s in Compare is not in Group info', variables = (diff_groups), code = "33704305")
        pass
    def set_resource(self):
        # cmp_num = len(self.option("cmp").prop['cmp_list'])
        # if cmp_num <= self.option("pool"):
        #     self.option("pool", cmp_num)
        # self._cpu = self.option("pool") + 1
        # self._memory = "{}G".format(self.option("pool")*30)
        self._cpu = 6
        self._memory = '50G'
    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "差异分析结果目录"]
            ])
        result_dir.add_regexp_rules([
            [r"*_vs_*.xls", "xls", "差异分析结果总表"],
            [r"*.DE.list", "xls", "差异基因列表"],
            [r"*summary.xls", "xls", "差异统计表"],
            ])
        super(DiffexpAgent, self).end()


class DiffexpTool(Tool):
    """
    Differential analysis based on edgeR, DEGseq, DEseq2.
    """
    def __init__(self, config):
        super(DiffexpTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'program/Python/bin/python'
        self.diff_toolbox = self.config.PACKAGE_DIR + "/tool_lab/diff_toolbox.py"
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
        self.rscript = '/program/R-3.3.1/bin/Rscript'
        self.r_volcano = self.config.PACKAGE_DIR + "/tool_lab/diff_volcano.r"


        group_dict = dict()
        with open(self.option('group').path, 'r') as group:
            for line in group.readlines():
                if line.startswith('#') or not line.strip():
                    continue
                else:
                    sample, group = line.strip().split('\t')
                    group_dict.setdefault(group, list())
                    group_dict[group].append(sample)
        print group_dict
        len_list = list()
        for g in group_dict:
            len_list.append(int(len(group_dict[g])))
            len_list.sort()

            # len_list.append(len(group_dict[g])).sort()
        if len_list[0] == len_list[-1] == 1:
            self.is_duplicate = 'False'
        else:
            self.is_duplicate = 'True'

    def diffexp(self):
        cmd = '{} {} '.format(self.python_path, self.diff_toolbox)
        cmd += '-count {} '.format(self.option("count").prop['path'])
        cmd += '-method {} '.format(self.option("method"))
        cmd += '-group {} '.format(self.option("group").prop['path'])
        cmd += '-cmp {} '.format(self.option("cmp").prop['path'])
        if self.option('norm'):
            cmd += '--normalized {} '.format(self.option('norm'))
        if self.option('model'):
            cmd += '-model {} '.format(self.option('model'))
        if self.option('test'):
            cmd += '-test {} '.format(self.option('test'))
        if self.option('voom'):
            cmd += '-voom {} '.format(self.option('voom'))
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
        if self.option('is_batch') == True and self.option('has_batch') == True:
            cmd += '--batch {} '.format(self.option('batch_matrix').prop['path'])
        if self.option('method').lower() == 'noiseq':
            cmd += '--is_duplicate {} '.format(self.is_duplicate)
            cmd += '--prob {} '.format(self.option('prob'))
        #cmd += '--tpm_filter_threshold {}'.format(self.option("tpm_filter_threshold"))
        print cmd
        cmd_name = 'diffexp'
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

    def volcano(self):
        static_graph_path = os.path.join(self.work_dir, 'static_graph')
        os.mkdir(static_graph_path)
        pvalue_padjust = self.option('pvalue_padjust')
        cmp_list = list()
        diff_files = glob.glob(os.path.join(self.work_dir, '*_vs_*.*.xls'))
        for each in diff_files:
            if each.endswith('.annot.xls'):
                continue
            if each.endswith('.normalize.xls'):
                continue
            if each.endswith('.sizeFactor.xls'):
                continue
            if each.endswith('.normFactor.xls'):
                continue
            diff_pd = pd.read_table(each, header=0, sep='\t')
            fname = os.path.basename(each)
            ctrl, test = re.match('(.*)_vs_(.*).{}.xls'.format(self.option('method').lower()), fname).groups()
            cmp_combine = ctrl + '|' + test
            cmp_list.append(cmp_combine)
            cmp_pd = pd.DataFrame([cmp_combine] * diff_pd.shape[0], columns=['compare'])
            if self.option('method').lower() in ['degseq', 'edger', 'deseq2', 'limma']:
                status_list, stat_cutoff, down_num, up_num = _get_volcano_status_cutoff(diff_pd, pvalue_padjust)
                volcano_pd = diff_pd.loc[:, ['seq_id', 'log2fc', pvalue_padjust, 'significant', 'regulate']]
                bool_ind = volcano_pd[pvalue_padjust] <= 0
                min_pvalue = min([x if x > 0 else '' for x in volcano_pd[pvalue_padjust].tolist()])
                volcano_pd.loc[bool_ind, pvalue_padjust] = min_pvalue
                volcano_pd[pvalue_padjust] = -np.log10(volcano_pd[pvalue_padjust])
                volcano_pd.dropna(inplace=True)
                regulate = list()
                for i in volcano_pd.index:
                    if volcano_pd['significant'].iloc[i] == 'no':
                        regulate.append(status_list[0])
                    if volcano_pd['significant'].iloc[i] == 'yes':
                        if volcano_pd['regulate'].iloc[i] == 'up':
                            regulate.append(status_list[2])
                        if volcano_pd['regulate'].iloc[i] == 'down':
                            regulate.append(status_list[1])
                volcano_pd['regulate'] = regulate
                volcano_pd.columns = ['seq_id', 'log2fc', 'log10pvalue', 'significant', 'regulate']
                bool_ind = volcano_pd['log10pvalue'] > stat_cutoff
                volcano_pd.loc[bool_ind, 'log10pvalue'] = stat_cutoff
            if self.option('method').lower() in ['noiseq']:
                status_list, down_num, up_num = _get_volcano_status_cutoff_noiseq(diff_pd)
                # sig_status[cmp_combine] = status_list
                volcano_pd = diff_pd.loc[:, ['seq_id', 'log2fc', 'D', 'significant', 'regulate']]
                volcano_pd.dropna(inplace=True)
                regulate = list()
                for i in volcano_pd.index:
                    if volcano_pd['significant'].iloc[i] == 'no':
                        regulate.append(status_list[0])
                    if volcano_pd['significant'].iloc[i] == 'yes':
                        if volcano_pd['regulate'].iloc[i] == 'up':
                            regulate.append(status_list[2])
                        if volcano_pd['regulate'].iloc[i] == 'down':
                            regulate.append(status_list[1])
                volcano_pd['regulate'] = regulate
                volcano_pd.columns = ['seq_id', 'log2fc', 'D', 'significant', 'regulate']
            volcano_pd = pd.concat([volcano_pd, cmp_pd], axis=1)
            volcano_pd_nosig = volcano_pd[volcano_pd['significant'] == 'no']
            # random select 10000 not sig diff genes for plotting
            if volcano_pd_nosig.shape[0] > 8000:
                volcano_pd_sig = volcano_pd[volcano_pd['significant'] == 'yes']
                volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.5)
                if volcano_pd_nosig.shape[0] > 12000:
                    volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.7)
                if volcano_pd_nosig.shape[0] > 12000:
                    volcano_pd_nosig = volcano_pd_nosig.sample(frac=0.7)
                volcano_pd = pd.concat([volcano_pd_sig, volcano_pd_nosig], axis=0)

            volcano_pd.to_csv(os.path.join(static_graph_path, '{}_vs_{}'.format(ctrl,test)), sep='\t',header=True,index=False)

    def static_volcano(self):
        cmd = '{} {} {} {}'.format(self.rscript, self.r_volcano, os.path.join(self.work_dir, 'static_graph'), self.option('method').lower())
        cmd += ' {}'.format(self.option('pvalue_padjust').lower())
        cmd_name = 'run_static_volcano'
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



    def set_output(self):
        volcano_path = os.path.join(self.work_dir, 'static_graph')
        files = glob.glob(os.path.join(volcano_path, '*.pdf'))
        for old_path in files:
            if old_path.endswith('.pdf'):
                file_name = os.path.basename(old_path)
                new_path = os.path.join(self.output_dir, file_name)
                os.link(old_path, new_path)
            else:
                continue
        if self.option("method").lower() in ["degseq", "edger", "deseq2", 'limma', 'svaseqlimma']:
            diff_files = glob.glob(self.option("output") + '/*_vs_*.*.xls')
            #diff_list = glob.glob(self.option("output") + '/*.DE.list')
            diff_summary = glob.glob(self.option("output") + '/*summary*.xls')
            diff_summary_pd = pd.read_table(diff_summary[0], header=0, sep="\t", index_col=0)
            diff_summary_pd_remove_no = diff_summary_pd[(True-diff_summary_pd['sum'].isin([0]))]
            diff_summary_pd_remove_no.to_csv(diff_summary[0], sep="\t")
            all_files = diff_files + diff_summary
            for each in all_files:
                fname = os.path.basename(each)
                if fname.endswith("normalize.xls"):
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
            #f_head = "seq_id\t \t \t"
            for diff_file_stat in diff_files_stat:
                ctrl, test = os.path.basename(diff_file_stat).split('.{}.xls'.format(self.option("method").lower()))[0].split('_vs_')
                fname = os.path.basename(diff_file_stat).split(".")[0]
                #f_head = f_head + fname + "\t*\t*\t*\t*\t*\t*\t"
                df_t = pd.read_table(diff_file_stat, index_col="seq_id")
                df_t.rename(columns={"log2fc": "{}_log2fc({}/{})".format(fname,test, ctrl), "fc": "{}_fc({}/{})".format(fname,test, ctrl),"pvalue":"{}_pvalue".format(fname),"padjust":"{}_padjust".format(fname),"significant":"{}_significant".format(fname),"regulate":"{}_regulate".format(fname)},
                              inplace=True)
                df_core = pd.DataFrame(df_t, columns=["{}_fc({}/{})".format(fname,test, ctrl),"{}_log2fc({}/{})".format(fname,test, ctrl),"{}_pvalue".format(fname), "{}_padjust".format(fname), "{}_significant".format(fname), "{}_regulate".format(fname)])
                df1 = pd.concat([df1, df_core], axis=1)
            df1.index.set_names("seq_id", inplace=True)
            df1.to_csv(self.output_dir + "/total_diff_stat.{}.xls".format(self.option("method").lower()), sep="\t")
            # if len(diff_files_stat)>1:
            #     with open(self.output_dir + "/total_diff_stat.{}.xls".format(self.option("method").lower()), 'r+') as f:
            #         content = f.read()
            #         f.seek(0, 0)
            #         f.write("seq_id" + content)
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
                df_t = pd.read_table(diff_file_stat, index_col="seq_id")
                df_t.rename(columns={"log2fc": "{}_log2fc({}/{})".format(fname, test, ctrl),
                                     "fc": "{}_fc({}/{})".format(fname, test, ctrl),
                                     "{}_mean".format(ctrl): "{}_{}_mean".format(fname,ctrl), "{}_mean".format(test): "{}_{}_mean".format(fname,test),
                                     "prob": "{}_prob".format(fname), "theta": "{}_theta".format(fname),
                                     "D": '{}_D'.format(fname),
                                     "significant": "{}_significant".format(fname),
                                     "regulate": "{}_regulate".format(fname)},
                            inplace=True)
                df_core = pd.DataFrame(df_t, columns=["{}_{}_mean".format(fname,ctrl), "{}_{}_mean".format(fname,test),
                                                      "{}_fc({}/{})".format(fname, test, ctrl),
                                                      "{}_log2fc({}/{})".format(fname, test, ctrl),
                                                      "{}_theta".format(fname),
                                                      '{}_D'.format(fname),
                                                      "{}_prob".format(fname),
                                                      "{}_significant".format(fname), "{}_regulate".format(fname)])
                df1 = pd.concat([df1, df_core], axis=1)
            df1.index.set_names("seq_id", inplace=True)
            df1.to_csv(self.output_dir + "/total_diff_stat.{}.xls".format(self.option("method").lower()), sep="\t")
            # if len(diff_files_stat)>1:
            #     with open(self.output_dir + "/total_diff_stat.{}.xls".format(self.option("method").lower()), 'r+') as f:
            #         content = f.read()
            #         f.seek(0, 0)
            #         f.write("seq_id" + content)
    def run(self):
        super(DiffexpTool, self).run()
        self.diffexp()
        self.volcano()
        self.static_volcano()
        self.set_output()
        self.end()

def _get_volcano_status_cutoff(diff_table, pvalue_padjust):
    sig_status = list()
    sig_mark = diff_table['significant']
    reg_list = diff_table['regulate']
    down_num = 0
    up_num = 0
    if 'no' in list(sig_mark):
        no_sig_num = sig_mark[sig_mark == "no"].shape[0]
        sig_status.append('nosig_' + str(no_sig_num))
    if 'yes' in list(sig_mark):
        reg_mark = reg_list[sig_mark == 'yes']
        if 'down' in list(reg_mark):
            down_num = reg_mark[reg_mark == 'down'].shape[0]
        sig_status.append('down_' + str(down_num))
        if 'up' in list(reg_mark):
            up_num = reg_mark[reg_mark == 'up'].shape[0]
        sig_status.append('up_' + str(up_num))

    sig_pvalues = diff_table[pvalue_padjust][diff_table['significant'] == "yes"]
    log10_sig_pvalues = -np.log10(sig_pvalues)
    log10_pvalue_list = sorted(list(log10_sig_pvalues[log10_sig_pvalues > 0]))

    if len(sig_pvalues) > 2000:
        log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list) * 0.85)]
    elif len(sig_pvalues) > 1000:
        log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list) * 0.90)]
    elif len(sig_pvalues) > 500:
        log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list) * 0.95)]
    elif len(sig_pvalues) > 250:
        log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list) * 0.99)]
    elif len(sig_pvalues) == 0:
        tmp = -np.log10(diff_table[pvalue_padjust])
        tmp_list = sorted(tmp[tmp > 0])
        if len(tmp_list) == 0:
            log10_pvalue_cutoff = 200
        else:
            log10_pvalue_cutoff = tmp_list[int(len(tmp_list) * 0.9)]
    else:
        # print(pvalue_padjust, diff_table, log10_pvalue_list)
        log10_pvalue_cutoff = log10_pvalue_list[int(len(log10_pvalue_list) * 0.8)]
    return sig_status, log10_pvalue_cutoff, down_num, up_num

def _get_volcano_status_cutoff_noiseq(diff_table):
    sig_status = list()
    sig_mark = diff_table['significant']
    reg_list = diff_table['regulate']
    down_num = 0
    up_num = 0
    if 'no' in list(sig_mark):
        no_sig_num = sig_mark[sig_mark == "no"].shape[0]
        sig_status.append('nosig_' + str(no_sig_num))
    if 'yes' in list(sig_mark):
        reg_mark = reg_list[sig_mark == 'yes']
        if 'down' in list(reg_mark):
            down_num = reg_mark[reg_mark == 'down'].shape[0]
        sig_status.append('down_' + str(down_num))
        if 'up' in list(reg_mark):
            up_num = reg_mark[reg_mark == 'up'].shape[0]
        sig_status.append('up_' + str(up_num))

    return sig_status, down_num, up_num

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
            "name": "tool_lab.diffexp",
            "instant": False,
            "options": dict(
                count='/mnt/ilustre/users/sanger-dev/workspace/20200113/Refrna_tsg_36819/Quant/output/gene.count.matrix',
                # exp='/mnt/ilustre/users/sanger-dev/workspace/20200113/Refrna_tsg_36819/Quant/output/gene.tpm.matrix',
                # method="edgeR",
                group="/mnt/ilustre/users/sanger-dev/workspace/20200113/Refrna_tsg_36819/remote_input/group_table/group.txt",
                cmp= "/mnt/ilustre/users/sanger-dev/workspace/20200113/Refrna_tsg_36819/remote_input/control_file/control.txt",
                pool=6,
                output=None,
                count_cutoff=4,
                dispersion=0.05,
                # padjust_way="BH",
                # pvalue=0.05,
                fc=2,
            )
        }
        # data['options']['method'] = 'edgeR'
        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.run()
        # #
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
        #
        # data['id'] += '3'
        # data['options']['method'] = 'limma'
        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.run()

        # data['id'] += '4'
        data['options']['method'] = 'NOIseq'
        data['options']['prob'] = '0.9'
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()


