# -*- coding: utf-8 -*-
# __author__ = 'zhaozhigang'
# @20210407

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil
import re
import subprocess
import itertools
import pandas as pd


class BugbaseContributionAgent(Agent):
    """
    多样性进行BugBase贡献度分析
    """
    def __init__(self, parent):
        super(BugbaseContributionAgent, self).__init__(parent)
        options = [
            {"name": "otu_table", "type": "infile", "format": "sequence.fasta"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "bugbase_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "species_level", "type": "int", "default": 9},
            {"name": "bugbase_id", "type": "string"},
            {"name": "otu_id", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "top", "type": "int", "default": 10},
            {"name": "method", "type": "string", "default": "average"},
        ]

        self.add_option(options)

    def step_start(self):
        self.step.bugbase_predict.start()
        self.step.update()

    def step_end(self):
        self.step.bugbase_predict.finish()
        self.step.update()

    def check_options(self):
        """
        进行参数二次检查
        :return:
        """
        if not self.option("otu_table").is_set:
            raise OptionError("请传入otu_table！")
        if not self.option("group_table").is_set:
            raise OptionError("请传入group_table！")
        if not self.option("bugbase_table").is_set:
            raise OptionError("请传入bugbase_table！")

    def set_resource(self):
        """
        设置资源
        :return:
        """
        self._cpu = 4
        self._memory = '20G'

    def end(self):
        """
        结束啦
        :return:
        """
        super(BugbaseContributionAgent, self).end()


class BugbaseContributionTool(Tool):
    def __init__(self, config):
        super(BugbaseContributionTool, self).__init__(config)
        """
        设置软件、脚本、数据库的路径
        """

    def run_profile_select(self):
        self.logger.info("start profile!")
        if self.option("method") == "none":
            options = {
                "origin_table": self.option("otu_table"),
                "samples":"all"
            }
        else:
            if self.option("method") == "average":
                method = 2
            else:
                method = 3
            options = {
                "group": self.option("group_table"),
                "origin_table": self.option("otu_table"),
                "group_method": method
            }
        self.gene_profile_select.set_options(options)
        self.gene_profile_select.on('end', self.run_contribute)
        self.gene_profile_select.run()


    def get_info(self):
        """
        根据代表序列获取otu与greengene数据库中的对应关系
        通过比对GreenGene--16s数据库，得到97%相似性的物种或者otu对应关系
        """
        self.logger.info("开始用picrust获取greengene数据库中的对应关系")
        if os.path.exists(os.path.join(self.work_dir, "Results")):
            shutil.rmtree(os.path.join(self.work_dir, "Results"))
        else:
            os.mkdir(os.path.join(self.work_dir, "Results"))
        result_path = os.path.join(self.work_dir, "Results")
        self.logger.info("设置picrust的结果文件完成")

        if not os.path.exists(os.path.join(self.work_dir,"otu_picking_params_97.txt")):
            self.logger.info("判断txt文件是否存在")
            f = open("otu_picking_params_97.txt", "wb")
            f.write("pick_otus:enable_rev_strand_match True\n")
            f.write("pick_otus:similarity 0.97")
            f.close()
        otu_params_path = os.path.join(self.work_dir, "otu_picking_params_97.txt")
        cmd = "{}pick_closed_reference_otus.py -i {} -o {} -f -p {} -r {}".format(self.python_scripts, self.option("otu_fasta").prop["path"], result_path, otu_params_path, self.Fundb)
        self.logger.info(cmd)
        command = self.add_command("picrust_predict", cmd).run()
        self.wait(command)
        if command.return_code in [0]:
            self.logger.info("picrust运行完成！")
        else:
            self.set_error("picrust预测失败！")

    def run(self):
        """
        运行
        """
        super(BugbaseContributionTool, self).run()
        self.logger.info("开始")
        self.run_profile_select()
        self.logger.info("运行tool结束")
        self.end()