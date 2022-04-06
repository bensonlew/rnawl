# -*- coding: utf-8 -*-
# __author__ = 'gudeqing, qinjincheng'

from biocluster.agent import Agent
import os
import glob
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
import pandas as pd


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
            dict(name="count", type="infile", format="small_rna.express_matrix"),
            dict(name="exp", type="infile", format="small_rna.express_matrix"),

            # group为样本分组信息文件，要求至少两列，第一列为样本名，其他列为组名。
            # 没有重复实验时，可以用样本名作为组名。
            dict(name="group", type="infile", format="small_rna.group_table"),
            # cmp为比较信息文件，仅两列，第一列为对照组名，第二列为实验组名，要求与group信息保持一致。
            dict(name="cmp", type="infile", format="small_rna.compare_table"),

            # count_cutoff 和 over_cutoff 是一起用来判定某个基因是否可以进入后续差异分析的的依据。
            # 对于一个基因，N个样本中至少有over_cutoff个样本的表达超过count_cutoff的才能进入后续分析。
            dict(name="count_cutoff", type="int", default=4),
            dict(name="over_cutoff", type="int", default=None),

            # 当用edgeR进行单样本和单样本比较时，需要输入经验值dispersion
            dict(name="dispersion", type="float", default=0.1),

            # pool为进程池的size. 用于差异分析的平行计算
            dict(name='pool', type='int', default=6),

            # pvalue 统计判断阈值
            dict(name="pvalue", type="float", default=0.05),
            dict(name="pvalue_padjust", type="string", default="padjust"),
            dict(name="fc", type="float", default=2),
            dict(name="padjust_way", type='string', default="BH"),

            #NOIseq独有的阈值
            dict(name='prob', type='float', default=0.8),

            # method: DESeq2, edgeR, DEGseq, limma, NOIseq
            dict(name="method", type="string", default="DESeq2"),

            # output
            dict(name="output", type="string", default=None),
            dict(name="diff_result", type="string"),
            dict(name="diff_list", type="string", ),
            dict(name="diff_summary", type="string", ),

        ]
        self.add_option(options)
        #当内存不足时，自动增加内存，每次增加30
        self._memory_increase_step = 30

    def check_options(self):
        if self.option("method").lower() not in ["degseq", "edger", "deseq2", 'noiseq', 'limma']:
            raise OptionError("method is incorrect")
        if self.option('method').lower() in ["degseq", "edger", "deseq2", 'limma']:
            pvalue, fc = float(self.option("pvalue")), float(self.option("fc"))
            dispersion = float(self.option("dispersion"))
            if not (0 <= pvalue <= 1):
                raise OptionError("pvalue cutoff must be in range [0-1]")
            if not (0 <= dispersion <= 1):
                raise OptionError("dispersion argument of edgeR must be in range [0-1]")
            if not fc >= 0:
                raise OptionError("fold change cutoff must be positive")
        groups = self.option("group").prop['group_dict'].keys()
        cmp_groups = [x for y in self.option("cmp").prop['cmp_list'] for x in y]
        diff_groups = set(cmp_groups) - set(groups)
        if diff_groups:
            raise OptionError('{} in compare is not in group info'.format(diff_groups))

    def set_resource(self):
        cmp_num = len(self.option("cmp").prop['cmp_list'])
        if cmp_num <= self.option("pool"):
            self.option("pool", cmp_num)
        self._cpu = self.option("pool") + 1
        self._memory = "{}G".format(10+self.option("pool")*6)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ['.', '', 'diffexp_tool_output_dir']
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
        self.python_path = 'miniconda2/bin/python'
        self.diff_toolbox = self.config.PACKAGE_DIR + "/small_rna_v2/diff_toolbox.py"
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
        # 判断是否有生物学重复,在method==noiseq时调用该参数
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

    def diffexp(self):
        cmd = '{} {} '.format(self.python_path, self.diff_toolbox)
        cmd += '-count {} '.format(self.option("count").prop['path'])
        cmd += '-exp {} '.format(self.option("exp").prop['path'])
        cmd += '--exp_type {} '.format('tpm')
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
        if self.option("over_cutoff") is not None:
            cmd += '--passed_number_cutoff {} '.format(self.option("over_cutoff"))
        if self.option('method').lower() == 'noiseq':
            cmd += '--is_duplicate {} '.format(self.is_duplicate)
            cmd += '--prob {} '.format(self.option('prob'))
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
                self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "32002706")
        else:
            self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "32002707")

    def set_output(self):
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
            diff_files_stat = glob.glob(os.path.join(self.output_dir, '*_vs_*.{}.xls'.format(self.option("method").lower())))
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
        #
        #
        # diff_files = glob.glob(self.option("output") + '/*_vs_*.*.xls')
        # # diff_list = glob.glob(self.option("output") + '/*.DE.list')
        # diff_summary = glob.glob(self.option("output") + '/*summary*.xls')
        # diff_summary_pd = pd.read_table(diff_summary[0], header=0, sep="\t", index_col=0)
        # diff_summary_pd_remove_no = diff_summary_pd[(True - diff_summary_pd['sum'].isin([0]))]
        # diff_summary_pd_remove_no.to_csv(diff_summary[0], sep="\t")
        # all_files = diff_files + diff_summary
        # for each in all_files:
        #     fname = os.path.basename(each)
        #     if fname.endswith("normalize.xls"):
        #         continue
        #     if fname.endswith("rawnormal.xls"):
        #         continue
        #     if fname.endswith("normFactor.xls"):
        #         continue
        #     if fname.endswith("sizeFactor.xls"):
        #         continue
        #     link = os.path.join(self.output_dir, fname)
        #     if os.path.exists(link):
        #         os.remove(link)
        #     os.link(each, link)

    def run(self):
        super(DiffexpTool, self).run()
        self.diffexp()
        self.set_output()
        self.end()

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):

        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import random

        exp_list = [
            {'count': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/known_miR_count.xls',
             'exp': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/known_miR_norm.xls'},
            {'count': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/novel_miR_count.xls',
             'exp': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/novel_miR_norm.xls'},
            {'count': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/all_miR_count.xls',
             'exp': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/exp/all_miR_norm.xls'},
        ]

        for exp_dict in exp_list:

            for method in ['edgeR', 'DESeq2', 'DEGseq']:

                data = {
                    'id': 'diffexp_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
                    'type': 'tool',
                    'name': 'small_rna.diffexp',
                    'instant': False,
                    'options': {
                        'count': exp_dict['count'],
                        'exp': exp_dict['exp'],
                        'method': method,
                        'group': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/group.txt',
                        'cmp': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/small_rna/control.txt',
                        'pool': 4,
                        'output': None,
                        'count_cutoff': 5,
                        'dispersion': 0.05,
                        'padjust_way': 'BH',
                        'pvalue': 0.05,
                        'fc': 2,
                    },
                }

                wsheet = Sheet(data=data)
                wf = SingleWorkflow(wsheet)
                wf.run()

if __name__ == '__main__':
    unittest.main()