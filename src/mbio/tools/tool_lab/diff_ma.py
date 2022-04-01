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
from mbio.packages.whole_transcriptome.utils import runcmd
__author__ = 'gdq'


class DiffMaAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(DiffMaAgent, self).__init__(parent)
        options = [
            # 输入表格 分三列,第一列为id,第二列为log2fc,第三列为pajust
            dict(name="project_type", type='string', default="custom"),
            dict(name="raw_file", type="infile", format="ref_rna_v2.common"),
            dict(name="pvalue", type="float", default=0.05),
            dict(name="fc", type="float", default=2),
            dict(name="x_axis_name", type='string', default="log10(TPM)"),
            dict(name="y_axis_name", type='string', default="log2(FC)"),
            dict(name="title_name", type='string', default="MA Plot"),
            dict(name="color", type='string', default="red_blue_grey")
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(DiffMaAgent, self).end()


class DiffMaTool(Tool):
    """
    Differential analysis based on edgeR, DEGseq, DEseq2.
    """
    def __init__(self, config):
        super(DiffMaTool, self).__init__(config)
        self.rscript = '/bioinfo/rna/miniconda2/bin/Rscript'
        self.r_volcano = self.config.PACKAGE_DIR + "/tool_lab/diff_ma2.r"


    def volcano(self):
        cmd = '{} {} '.format(self.rscript, self.r_volcano)
        cmd += '-i {} '.format(os.path.join(self.work_dir,"plot.txt"))
        cmd += '-c {} '.format(os.path.join(self.work_dir,"plot_config.txt"))
        print cmd
        cmd_name = 'diff_ma'
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
        if self.option("project_type") == "custom" :
            if self.option("raw_file").prop["path"].endswith("txt"):
                raw_df = pd.read_table(self.option("raw_file").prop["path"])
            elif self.option("raw_file").prop["path"].endswith("xls"):
                try:
                    raw_df = pd.read_table(self.option("raw_file").prop["path"])
                except:
                    raw_df = pd.read_excel(self.option("raw_file").prop["path"])
            else:
                self.set_error('输入文件必须为txt或者xls格式之一')
            raw_df.columns = ["seq_id", "group1_tpm", "group2_tpm", "log2fc", "padjust"]
            raw_df["log10fpkm"] = raw_df.loc[:, ["group1_tpm", "group2_tpm"]].mean(axis=1)
            raw_df = raw_df[raw_df["log10fpkm"] != 0]
            raw_df['log10fpkm'] = (raw_df["log10fpkm"]).apply(np.log10)
            raw_df["regulate"] = raw_df.apply(lambda x: "up" if (x["log2fc"] > math.log(self.option("fc"),2) and x["padjust"] < self.option("pvalue")) else "down" if (
                        x["log2fc"] < -math.log(self.option("fc"),2) and x["padjust"] < self.option("pvalue")) else "normal", axis=1)
            raw_df = raw_df[["seq_id", "log10fpkm", "log2fc", "regulate"]]
            raw_df.to_csv(os.path.join(self.work_dir,"plot.txt"),sep="\t",index=False)
        elif self.option("project_type") == "ref_rna_v2" :
            try:
                diff_pd = pd.read_table(self.option("raw_file").prop["path"], header=0, sep='\t')
                columns = diff_pd.columns
                fc_ind = list(columns).index('log2fc')
                need_cols = [columns[0]]
                need_cols = need_cols + [columns[fc_ind - 3], columns[fc_ind - 2]]
                need_cols += ['log2fc', "padjust"]
                need_df = diff_pd.loc[:, need_cols]
                need_df["log10fpkm"] = need_df.iloc[:, [1, 2]].mean(axis=1)
                need_df = need_df[need_df["log10fpkm"] != 0]
                need_df['log10fpkm'] = (need_df["log10fpkm"]).apply(np.log10)
                need_df["regulate"] = need_df.apply(lambda x: "up" if (
                            x["log2fc"] > math.log(self.option("fc"), 2) and x["padjust"] < self.option(
                        "pvalue")) else "down" if (
                        x["log2fc"] < -math.log(self.option("fc"), 2) and x["padjust"] < self.option(
                    "pvalue")) else "normal", axis=1)
                need_df = need_df[[columns[0], "log10fpkm", "log2fc", "regulate"]]
                need_df.to_csv(os.path.join(self.work_dir, "plot.txt"), sep="\t", index=False)
            except:
                diff_pd = pd.read_table(self.option("raw_file").prop["path"], header=0, sep='\t')
                columns = diff_pd.columns
                for col in columns:
                    if "log2fc" in col:
                        select_col = col
                fc_ind = list(columns).index(select_col)
                need_cols = [columns[0]]
                tpm_list = []
                for col in columns:
                    if col.endswith("tpm") or col.endswith("fpkm"):
                        tpm_list.append(col)
                need_cols = need_cols + [tpm_list[-2], tpm_list[-1]]
                need_cols += [select_col, "padjust"]
                need_df = diff_pd.loc[:, need_cols]
                need_df["log10fpkm"] = need_df.iloc[:, [1, 2]].mean(axis=1)
                need_df = need_df[need_df["log10fpkm"] != 0]
                need_df['log10fpkm'] = (need_df["log10fpkm"]).apply(np.log10)
                need_df["regulate"] = need_df.apply(lambda x: "up" if (
                        x[select_col] > math.log(self.option("fc"), 2) and x["padjust"] < self.option(
                    "pvalue")) else "down" if (
                        x[select_col] < -math.log(self.option("fc"), 2) and x["padjust"] < self.option(
                    "pvalue")) else "normal", axis=1)
                need_df = need_df[[columns[0], "log10fpkm", select_col, "regulate"]]
                need_df.columns = [columns[0], "log10fpkm", "log2fc", "regulate"]
                need_df.to_csv(os.path.join(self.work_dir, "plot.txt"), sep="\t", index=False)




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
        config_dict["x"] = self.option("x_axis_name")
        config_dict["y"] = self.option("y_axis_name")
        config_dict["t"] = self.option("title_name")
        config_df = pd.DataFrame([config_dict])
        config_df.to_csv(os.path.join(self.work_dir,"plot_config.txt"),sep="\t",index=False)



    def set_output(self):
        for i in ["diff_ma.pdf","diff_ma.png",'diff_ma.svg']:
            if os.path.exists(os.path.join(self.output_dir,i)):
                os.remove(os.path.join(self.output_dir,i))
            os.link(os.path.join(self.work_dir,i),os.path.join(self.output_dir,i))


    def run(self):
        super(DiffMaTool, self).run()
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
            "id": "DiffMa" + str(random.randint(1, 10000)),
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


