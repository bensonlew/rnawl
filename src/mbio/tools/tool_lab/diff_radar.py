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
from collections import OrderedDict
import numpy as np
from mbio.packages.whole_transcriptome.utils import runcmd
__author__ = 'gdq'


class DiffRadarAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(DiffRadarAgent, self).__init__(parent)
        options = [
            # 输入表格 分三列,第一列为id,第二列为log2fc,第三列为pajust
            dict(name="project_type", type='string', default="custom"),
            dict(name="raw_file", type="infile", format="ref_rna_v2.common"),
            dict(name="group_table", type="infile", format="ref_rna_v2.common"),
            dict(name="target_ids", type='string', default=""),
        ]
        self.add_option(options)

    def check_options(self):
        if  self.option("project_type") == "custom" and not self.option("group_table").is_set:
            self.set_error("当客户自主上传文件时,需提供group表格")
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(DiffRadarAgent, self).end()


class DiffRadarTool(Tool):
    """
    Differential analysis based on edgeR, DEGseq, DEseq2.
    """
    def __init__(self, config):
        super(DiffRadarTool, self).__init__(config)
        self.rscript = '/bioinfo/rna/miniconda2/bin/Rscript'
        self.r_volcano = self.config.PACKAGE_DIR + "/tool_lab/diff_ma2.r"


    def prepare_plot_file(self):
        if self.option("project_type") == "custom" :
            # if self.option("raw_file").prop["path"].endswith("txt"):
            #     raw_df = pd.read_table(self.option("raw_file").prop["path"])
            # elif self.option("raw_file").prop["path"].endswith("xls"):
            #     try:
            #         raw_df = pd.read_table(self.option("raw_file").prop["path"])
            #     except:
            #         raw_df = pd.read_excel(self.option("raw_file").prop["path"])
            # else:
            #     self.set_error('输入文件必须为txt或者xls格式之一')
            raw_df =  pd.read_table(self.option("raw_file").prop["path"],index_col= 0)
            if self.option("target_ids") :
                target_ids = self.option("target_ids").split(",")
                for id in target_ids:
                    if not id in raw_df.index:
                        self.set_error("输入id{}在差异文件中无法找到".format(id))
                plot_df = raw_df.loc[target_ids,]
            else:
                if raw_df.shape[0] > 30:
                    plot_df = raw_df.head(30)
                else:
                    plot_df = raw_df
            lo2fc_col = list(plot_df.columns)[-1]
            group_dict = OrderedDict()
            with open(self.option("group_table").prop["path"]) as r:
                # r.readline()
                for line in r.readlines():
                    line = line.strip().split("\t")
                    group_dict.setdefault(line[1], list())
                    if line[0] not in group_dict[line[1]]:
                        group_dict[line[1]].append(line[0])
            control_group = list(group_dict)[0]
            compare_group = list(group_dict)[1]
            plot_df[compare_group] = plot_df.loc[:, [i for i in group_dict[compare_group]]].mean(axis=1)
            plot_df[control_group] = plot_df.loc[:, [i  for i in group_dict[control_group]]].mean(axis=1)
            target_col = [control_group, compare_group, lo2fc_col]
            final_plot_df = plot_df[target_col]
            final_plot_df.to_csv(os.path.join(self.output_dir,"plot.txt"),sep="\t")
        elif self.option("project_type") == "ref_rna_v2" :
            try:
                diff_pd = pd.read_table(self.option("raw_file").prop["path"], header=0, sep='\t')
                columns = diff_pd.columns
                fc_ind = list(columns).index('log2fc')
                need_cols = [columns[0]]
                need_cols = need_cols + [columns[fc_ind - 3], columns[fc_ind - 2]]
                need_cols += ['log2fc']
                if self.option("target_ids"):
                    target_ids = self.option("target_ids").split(",")
                    for id in target_ids:
                        if not id in diff_pd[columns[0]].values:
                            self.set_error("输入id{}在差异文件中无法找到".format(id))
                    plot_df = diff_pd[diff_pd[columns[0]].isin(target_ids)]
                else:
                    if diff_pd.shape[0] > 30:
                        plot_df = diff_pd.head(30)

                need_df = plot_df.loc[:, need_cols]
                need_df.columns = [i.replace("_tpm", "") for i in need_df.columns]
                need_df.columns = [i.replace("_fpkm", "") for i in need_df.columns]
                # need_df = need_df[[columns[0], "log10fpkm", "log2fc", "regulate"]]
                need_df.to_csv(os.path.join(self.output_dir, "plot.txt"), sep="\t", index=False)
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
                need_cols += [select_col]
                if self.option("target_ids"):
                    target_ids = self.option("target_ids").split(",")
                    for id in target_ids:
                        if not id in diff_pd[columns[0]].values:
                            self.set_error("输入id{}在差异文件中无法找到".format(id))
                    plot_df = diff_pd[diff_pd[columns[0]].isin(target_ids)]
                else:
                    if diff_pd.shape[0] > 50:
                        plot_df = diff_pd.head(50)
                need_df = plot_df.loc[:, need_cols]
                need_df.columns = [i.replace("_tpm", "") for i in need_df.columns]
                need_df.columns =  [i.replace("_fpkm", "") for i in need_df.columns]
                need_df.to_csv(os.path.join(self.output_dir, "plot.txt"), sep="\t", index=False)
        else:
            if os.path.exists(os.path.join(self.output_dir, "plot.txt")):
                os.remove(os.path.join(self.output_dir, "plot.txt"))
            os.link(self.option("raw_file").prop["path"],os.path.join(self.output_dir, "plot.txt"))



    def run(self):
        super(DiffRadarTool, self).run()
        self.prepare_plot_file()
        # self.prepare_plot_config()
        # self.volcano()
        # self.set_output()
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
            "id": "DiffRadar" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.diff_radar",
            "instant": False,
            "options": dict(
                raw_file='s3://refrnav2/files/m_188/188_5ffbaeead3e00/mbs6_v5o6eq0967dj319dmpsmrq/workflow_results/07DiffExpress_G/HFL_vs_HGL.edger.annot.xls',
                # exp='/mnt/ilustre/users/sanger-dev/workspace/20200113/Refrna_tsg_36819/Quant/output/gene.tpm.matrix',
                # method="edgeR",
                project_type = "ref_rna_v2",
            )
        }

        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()


