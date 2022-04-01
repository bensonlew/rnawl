# -*- coding: utf-8 -*-
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import unittest
import re
import math
import numpy as np
import shutil
from mbio.packages.whole_transcriptome.utils import runcmd
__author__ = 'gdq'


class MegenaAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(MegenaAgent, self).__init__(parent)
        options = [
            # 输入表格 分三列,第一列为id,第二列为log2fc,第三列为pajust
            dict(name="project_type", type='string', default="custom"),
            dict(name="raw_file", type="infile", format="ref_rna_v2.common"),
            dict(name="corr_method", type="string", default="pearson"),
            dict(name="fdr_cutoff", type="float", default=0.05),
            dict(name="module_pval", type='float', default=0.05),
            dict(name="min_size", type='int', default=10)
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(MegenaAgent, self).end()


class MegenaTool(Tool):
    """
    """
    def __init__(self, config):
        super(MegenaTool, self).__init__(config)
        self.rscript = '/bioinfo/tool_lab/MEGENA/miniconda3/bin/Rscript'
        self.r_megena = self.config.PACKAGE_DIR + "/tool_lab/megena.r"


    def run_megena(self):
        if os.path.exists(os.path.join(self.work_dir,"tmp_dir")):
            shutil.rmtree(os.path.join(self.work_dir,"tmp_dir"))
        os.makedirs(os.path.join(self.work_dir,"tmp_dir"))
        cmd = '{} {} '.format(self.rscript, self.r_megena)
        cmd += '-e {} '.format(os.path.join(self.work_dir,"final_express.txt"))
        cmd += '-c {} '.format(os.path.join(self.work_dir,"plot_config.txt"))
        cmd += '-t {} '.format(os.path.join(self.work_dir,"tmp_dir"))
        print cmd
        cmd_name = 'megena'
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
                self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd))
        else:
            self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd))

    def prepare_plot_file(self):
        if self.option("project_type") == "custom" :
            final_express_matrix = os.path.join(self.work_dir,"final_express.txt")
            with open(self.option("raw_file").prop["path"],"r") as r,open(final_express_matrix,"w") as w:
                header = r.readline()
                new_header = "\t".join(header.strip().split("\t")[1:])+"\n"
                w.write(new_header)
                for line in r:
                    w.write(line)




    def prepare_plot_config(self):
        config_dict ={}
        config_dict["method"] = self.option("corr_method")
        config_dict["fdr_cutoff"] = self.option("fdr_cutoff")
        config_dict["module_pval"] = self.option("module_pval")
        config_dict["min_size"] = self.option("min_size")
        config_df = pd.DataFrame([config_dict])
        config_df.to_csv(os.path.join(self.work_dir,"plot_config.txt"),sep="\t",index=False)



    def set_output(self):
        # for i in ["diff_ma.pdf","diff_ma.png",'diff_ma.svg']:
        #     if os.path.exists(os.path.join(self.output_dir,i)):
        #         os.remove(os.path.join(self.output_dir,i))
        #     os.link(os.path.join(self.work_dir,i),os.path.join(self.output_dir,i))
        pass


    def run(self):
        super(MegenaTool, self).run()
        self.prepare_plot_file()
        self.prepare_plot_config()
        self.run_megena()
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
            "id": "Megena" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.diff_ma_new",
            "instant": False,
            "options": dict(
                raw_file='s3://refrnav2/files/m_188/188_5ffbaeead3e00/mbs6_v5o6eq0967dj319dmpsmrq/workflow_results/07DiffExpress_G/HFL_vs_HGL.edger.annot.xls',
                # exp='/mnt/ilustre/users/sanger-dev/workspace/20200113/Refrna_tsg_36819/Quant/output/gene.tpm.matrix',
                # method="edgeR",
                project_type = "ref_rna_v2",
                pvalue=0.05,
                fc= 2,
                x_axis_name="log10(TPM)",
                y_axis_name= "log2(FC)",
                title_name="MA Plot",
                color = "ref_blue_grey"
            )
        }

        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()


