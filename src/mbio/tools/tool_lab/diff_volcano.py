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


class DiffVolcanoAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(DiffVolcanoAgent, self).__init__(parent)
        options = [
            # 输入表格 分三列,第一列为id,第二列为log2fc,第三列为pajust
            dict(name="raw_file", type="infile", format="ref_rna_v2.common"),
            dict(name="pvalue", type="float", default=0.05),
            dict(name="fc", type="float", default=2),
            dict(name="x_axis_name", type='string', default="log2(FC)"),
            dict(name="y_axis_name", type='string', default="-log10(Padjust)"),
            dict(name="title_name", type='string', default="Volcano Plot"),
            dict(name="color", type='string', default="red_blue_grey"),
            dict(name="showname_num", type='string', default=""),
            dict(name="showname_str", type='string', default="")
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(DiffVolcanoAgent, self).end()


class DiffVolcanoTool(Tool):
    """
    Differential analysis based on edgeR, DEGseq, DEseq2.
    """
    def __init__(self, config):
        super(DiffVolcanoTool, self).__init__(config)
        self.rscript = '/bioinfo/rna/miniconda2/bin/Rscript'
        self.r_volcano = self.config.PACKAGE_DIR + "/tool_lab/diff_volcano_new2.r"


    def volcano(self):
        cmd = '{} {} '.format(self.rscript, self.r_volcano)
        cmd += '-i {} '.format(os.path.join(self.work_dir,"plot.txt"))
        cmd += '-c {} '.format(os.path.join(self.work_dir,"plot_config.txt"))
        if self.option("showname_num"):
            cmd += '-n {} '.format(self.option("showname_num"))
        if self.option("showname_str"):
            cmd += '-s {} '.format(self.option("showname_str"))
        print cmd
        cmd_name = 'diff_volcano'
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

    def prepare_plot_file(self):
        if self.option("raw_file").prop["path"].endswith("txt"):
            raw_df = pd.read_table(self.option("raw_file").prop["path"])
            # raw_df.iloc[:, 1][np.isinf(raw_df.iloc[:, 1])] = 0
        elif self.option("raw_file").prop["path"].endswith("xls"):
            try:
                raw_df = pd.read_table(self.option("raw_file").prop["path"])
            except:
                raw_df = pd.read_excel(self.option("raw_file").prop["path"])
        else:
            self.set_error('输入文件必须为txt或者xls格式之一')
        raw_df.iloc[:, 1][np.isinf(raw_df.iloc[:, 1])] = 0
        raw_df.columns = ["seq_id", "log2fc", "padjust"]
        raw_df.to_csv(os.path.join(self.work_dir,"plot.txt"),sep="\t",index=False)

    def prepare_plot_config(self):
        method = 1
        if self.option("color") == "red_blue_grey":
            method = 1
        elif self.option("color") == "red_green_grey":
            method = 2
        elif self.option("color") == "red_blue_black":
            method = 3
        elif self.option("color") == "red_green_black":
            method = 4
        config_dict ={}
        config_dict["m"] = method
        config_dict["l"] = self.option("fc")
        config_dict["p"] = self.option("pvalue")
        config_dict["x"] = self.option("x_axis_name")
        config_dict["y"] = self.option("y_axis_name")
        config_dict["t"] = self.option("title_name")
        config_df = pd.DataFrame([config_dict])
        config_df.to_csv(os.path.join(self.work_dir,"plot_config.txt"),sep="\t",index=False)



    def set_output(self):
        for i in ["volcano.pdf","volcano.png",'volcano.svg']:
            if os.path.exists(os.path.join(self.output_dir,i)):
                os.remove(os.path.join(self.output_dir,i))
            os.link(os.path.join(self.work_dir,i),os.path.join(self.output_dir,i))


    def run(self):
        super(DiffVolcanoTool, self).run()
        self.prepare_plot_file()
        self.prepare_plot_config()
        self.volcano()
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
            "id": "DiffVolcano" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.diff_volcano",
            "instant": False,
            "options": dict(
                # raw_file='/mnt/ilustre/users/sanger-dev/workspace/20210322/DiffexpBatch_jgo0_k3guobgpi1lrvksgs4pu61_6989_3353/DiffexpBatch/output/voca.xls',
                # exp='/mnt/ilustre/users/sanger-dev/workspace/20200113/Refrna_tsg_36819/Quant/output/gene.tpm.matrix',
                raw_file = "s3://toollab/files/test_diff_ma/test_diff_ma/diff_ma.pdf",
                # method="edgeR",
                pvalue=0.05,
                fc= 2,
                x_axis_name="log2(FC)",
                y_axis_name= "-log10(Padjust)",
                title_name="Volcano Plot",
                color = "ref_blue_grey"
            )
        }

        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()


