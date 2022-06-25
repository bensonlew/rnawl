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
import json
__author__ = 'gdq'


class PathwayNetworkAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(PathwayNetworkAgent, self).__init__(parent)
        options = [
            {"name": "raw_file", "type": "infile", "format": "ref_rna_v2.common"},  # fasta文件
            # {"name": "raw_file_info", "type": "string"},  # fasta文件
            {"name": "project_type", "type": "string", "default": None},  # 产品类型 [refrna,custom]等
            {"name": "pvalue_pajust", "type": "string", "default": "padjust"},
            {"name": "anno_type", "type": "string", "default": None},  # 注释类型 [kegg,go]
            {"name": "target_ids", "type": "string", "default": None},  # 目标id,以英文逗号分割
            {"name": "sup_connect", "type": "bool", "default": False},  # 目标id,以英文逗号分割
            {"name": "sup_all", "type": "bool", "default": False},  # 目标id,以英文逗号分割
        ]
        self.add_option(options)



    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(PathwayNetworkAgent, self).end()


class PathwayNetworkTool(Tool):
    """
    """
    def __init__(self, config):
        super(PathwayNetworkTool, self).__init__(config)
        if self.option("anno_type").lower() == "kegg":
            self.anno_dir = os.path.join(self.config.SOFTWARE_DIR,"database","Tool_lab","pathway_network","kegg")
        else:
            self.anno_dir =  os.path.join(self.config.SOFTWARE_DIR,"database","Tool_lab","pathway_network","go")
        self.python_path = "/miniconda2/bin/python"
        self.network_file = os.path.join(self.config.PACKAGE_DIR,"tool_lab","network_pathway_ana.py")
        python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/'
        self.set_environ(PATH=python_path)




    def prepare_plot_file(self):
        if self.option("project_type") == "custom" :
            raw_df =  pd.read_table(self.option("raw_file").prop["path"],index_col= 0)
            if self.option("target_ids") :
                target_ids = self.option("target_ids").split(",")
                for id in target_ids:
                    if not id in raw_df.index:
                        self.set_error("输入id{}在差异文件中无法找到".format(id))
                plot_df = raw_df.loc[target_ids,]
            else:
                if raw_df.shape[0] > 20:
                    plot_df = raw_df.head(20)
                else:
                    plot_df = raw_df
            analysis_df = plot_df
            analysis_df.columns = ["gene_num","padjust"]
            analysis_df.to_csv(os.path.join(self.work_dir, "prepare.txt"), sep="\t")
        elif self.option("project_type") == "ref_rna_v2" :
            if self.option("anno_type").lower() == "kegg":
                try:
                    raw_df = pd.read_table(self.option("raw_file").prop["path"])
                    if self.option("pvalue_pajust").lower() == "padjust":
                        b = raw_df.iloc[:, [3, 0, 7]]
                    else:
                        b = raw_df.iloc[:, [3, 0, 6]]
                    uni_df = b.set_index(b.columns[0])
                    if self.option("target_ids"):
                        target_ids = self.option("target_ids").split(",")
                        for id in target_ids:
                            if not id in uni_df.index:
                                self.logger.set_error("输入id{}在差异文件中无法找到".format(id))
                        plot_df = uni_df.loc[target_ids,]
                    else:
                        if uni_df.shape[0] > 20:
                            plot_df = uni_df.head(20)
                        else:
                            plot_df = uni_df
                    analysis_df = plot_df
                    analysis_df.columns = ["gene_num", "padjust"]
                    analysis_df.to_csv(os.path.join(self.work_dir, "prepare.txt"), sep="\t")
                except:
                    self.set_error("文件预处理失败")
            if self.option("anno_type").lower() == "go":
                a = pd.read_table(self.option("raw_file").prop["path"])
                if self.option("pvalue_pajust").lower() == "padjust":
                    b = a.iloc[:, [0, 11, 7]]
                else:
                    b = a.iloc[:, [0, 11, 6]]
                uni_df = b.set_index(b.columns[0])
                if self.option("target_ids"):
                    target_ids = self.option("target_ids").split(",")
                    for id in target_ids:
                        if not id in uni_df.index:
                            self.logger.set_error("输入id{}在差异文件中无法找到".format(id))
                    plot_df = uni_df.loc[target_ids,]
                else:
                    if uni_df.shape[0] > 20:
                        plot_df = uni_df.head(20)
                    else:
                        plot_df = uni_df
                analysis_df = plot_df
                analysis_df.columns = ["gene_num", "padjust"]
                analysis_df.to_csv(os.path.join(self.work_dir, "prepare.txt"), sep="\t")
        prepare_file = os.path.join(self.work_dir, "prepare.txt")
        anno_type = self.option("anno_type")
        sup_connect = self.option("sup_connect")
        sup_all = self.option("sup_all")
        json_dict = {
            "prepare_file": prepare_file,
            "anno_type" : anno_type ,
            "sup_connect" : sup_connect,
            "sup_all" : sup_all,
            "anno_dir" : self.anno_dir
        }
        with open(self.work_dir + "/pathway_network.json", 'w') as json_f:
            json.dump(json_dict, json_f, sort_keys=True, indent=4)


    def run_pathway_network(self):
        cmd = "{} {} ".format(self.python_path , self.network_file)
        # cmd = "{} {} ".format("python", self.network_file)
        cmd += "-json_file {} ".format(os.path.join(self.work_dir,"pathway_network.json"))
        cmd += "-output_dir {}".format(self.output_dir)
        cmd_name = 'pathway_network_analysis'
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
                self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd))
        else:
            self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd))




    def run(self):
        super(PathwayNetworkTool, self).run()
        self.prepare_plot_file()
        self.run_pathway_network()
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
            "id": "PathwayNetwork" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.pathway_network",
            "instant": False,
            "options": dict(
                raw_file='/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab_batch2/path_network/test_pathway.txt',
                # exp='/mnt/ilustre/users/sanger-dev/workspace/20200113/Refrna_tsg_36819/Quant/output/gene.tpm.matrix',
                anno_type="KEGG",
                project_type = "custom",
                sup_connect = True,
                sup_all =True
            )
        }

        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()


