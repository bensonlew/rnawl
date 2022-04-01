# -*- coding: utf-8 -*-
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
__author__ = 'gdq'


class ExpresscorAgent(Agent):
    def __init__(self, parent):
        super(ExpresscorAgent, self).__init__(parent)
        options = [
            dict(name="exp", type="infile", format="ref_rna_v2.express_matrix"),
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

            # method: DESeq2, edgeR, DEGseq
            dict(name="method", type="string", default="DESeq2"),

            # output
            dict(name="output", type="string", default=None),
            dict(name="diff_result", type="string"),
            dict(name="diff_list", type="string", ),
            dict(name="diff_summary", type="string", ),
        ]
        self.add_option(options)

    def check_options(self):
        if self.option("method").lower() not in ["degseq", "edger", "deseq2"]:
            raise OptionError("Method is incorrect")
        pvalue, fc = float(self.option("pvalue")), float(self.option("fc"))
        dispersion = float(self.option("dispersion"))
        if not (0 <= pvalue <= 1):
            raise OptionError("pvalue cutoff must be in range [0-1]")
        if not (0 <= dispersion <= 1):
            raise OptionError("Dispersion argument of edgeR must be in range [0-1]")
        if not fc >= 0:
            raise OptionError("fold change cutoff must be positive")
        groups = self.option("group").prop['group_dict'].keys()
        cmp_groups = [x for y in self.option("cmp").prop['cmp_list'] for x in y]
        diff_groups = set(cmp_groups) - set(groups)
        if diff_groups:
            raise OptionError('Groups:{} in Compare is not in Group info'.format(diff_groups))

    def set_resource(self):
        cmp_num = len(self.option("cmp").prop['cmp_list'])
        if cmp_num <= self.option("pool"):
            self.option("pool", cmp_num)
        self._cpu = self.option("pool") + 1
        self._memory = "{}G".format(self.option("pool")*30)

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
        super(ExpresscorAgent, self).end()


class ExpresscorTool(Tool):
    """
    Differential analysis based on edgeR, DEGseq, DEseq2.
    """
    def __init__(self, config):
        super(ExpresscorTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'program/Python/bin/python'
        self.diff_toolbox = self.config.PACKAGE_DIR + "/denovo_rna_v2/diff_toolbox.py"
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

    def Expresscor(self):
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
        cmd_name = 'Expresscor'
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
                self.set_error("运行{}>>>{}出错".format(cmd_name, cmd))
        else:
            self.set_error("运行{}>>>{}出错".format(cmd_name, cmd))

    def set_output(self):
        diff_files = glob.glob(self.option("output") + '/*_vs_*.xls')
        diff_list = glob.glob(self.option("output") + '/*.DE.list')
        diff_summary = glob.glob(self.option("output") + '/*summary.xls')
        all_files = diff_files + diff_list + diff_summary
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(ExpresscorTool, self).run()
        self.Expresscor()
        self.set_output()
        self.end()


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
            "id": "Expresscor" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.Expresscor",
            "instant": False,
            "options": dict(
                count=self.test_dir + '/transcript.count.matrix',
                exp=self.test_dir + '/transcript.tpm.matrix',
                method="edgeR",
                group=self.test_dir + "/default_group.txt",
                cmp=self.test_dir + "/control_file.txt",
                pool=6,
                output=None,
                count_cutoff=5,
                dispersion=0.05,
                padjust_way="BH",
                pvalue=0.05,
                fc=2,
            )
        }
        data['options']['method'] = 'edgeR'
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()
        #
        data['id'] += '1'
        data['options']['method'] = 'DESeq2'
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()
        #
        data['id'] += '2'
        data['options']['method'] = 'DEGseq'
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()


