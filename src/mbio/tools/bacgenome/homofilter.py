#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os, re
from mbio.packages.bacgenome.common import sum_stat
import unittest
import pandas as pd

class HomofilterAgent(Agent):
    """
    核心基因和特有基因
    version 1.0
    author: ysh
    last_modify: 2019.04.16
    """

    def __init__(self, parent):
        super(HomofilterAgent, self).__init__(parent)
        options = [
            {"name": "anno_dir", "type": "infile", "format": "annotation.mg_anno_dir"},
            {"name": "stat_file", "type": "infile", "format": "sequence.profile_table"},
            {'name': 'filter', 'type': 'string'},  # 过滤的样本
            {'name': 'select', 'type': 'string'},  # 筛选的样本
            {"name": "result", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('anno_dir').is_set:
            raise OptionError("必须设置anno文件目录！")
        if not self.option('stat_file').is_set:
            raise OptionError("必须输入同源基因stat文件！")
        if not self.option("filter") and not self.option("select"):
            raise OptionError("不能同时没有过滤样本和筛选样本！")
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        super(HomofilterAgent, self).end()


class HomofilterTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(HomofilterTool, self).__init__(config)
        self.python_path = 'miniconda2/bin/python'
        self.script = self.config.PACKAGE_DIR + "/bacgenome/homofilter.py"

    def run_homofilter(self):
        stat_file = self.option("stat_file").prop["path"]
        outfile = os.path.join(self.output_dir, "Homofilter_genes.xls")
        cmd = '{} {} -i {} -o {} -anno {}'.format(self.python_path, self.script, stat_file, outfile, self.anno)
        if self.option("select"):
            cmd += ' -s {} '.format(self.option("select"))
        if self.option("filter"):
            cmd += ' -f {} '.format(self.option("filter"))
        self.logger.info(cmd)
        command = self.add_command("homofilter", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_homofilter运行完成")
        else:
            self.set_error("run_homofilter运行出错!")

    def creat_anno(self):
        self.anno = os.path.join(self.work_dir, "anno.xls")
        anno_dir = self.option("anno_dir").prop["path"]
        all_files = os.listdir(anno_dir)
        table_list = []
        for each in all_files:
            file_path = os.path.join(anno_dir,each)
            sample = each.split(".xls")[0]
            table = pd.read_table(file_path, sep="\t", header=0)
            table["Sample"] = sample
            table["SampleGene"] = table["Sample"].str.cat(table["Gene ID"], sep="|")
            del table["Sample"]
            table_list.append(table)
        final_anno = pd.concat(table_list, axis=0)
        final_anno.to_csv(self.anno, sep="\t",index=False)

    def set_output(self):
        homofilter_file = os.path.join(self.output_dir, "Homofilter_genes.xls")
        self.option("result").set_path(homofilter_file)
        self.end()

    def run(self):
        """
        运行
        """
        super(HomofilterTool, self).run()
        self.creat_anno()
        self.run_homofilter()
        self.set_output()
