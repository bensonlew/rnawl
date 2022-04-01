# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os
import linecache
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class LefseSortSampleMgAgent(Agent):
    """
    传入一个group表，以及是否进行样本合并的参数生成一张丰度表并对并依照group表OTU表进行筛选合并
    """
    def __init__(self, parent):
        super(LefseSortSampleMgAgent, self).__init__(parent)
        options = [
            {"name": "in_otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入丰度文件
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 输入的group表
            {"name": "out_otu_table", "type": "outfile", "format": "meta.otu.otu_table"}  # 输出的结果OTU表
        ]
        self.add_option(options)
        self.step.add_steps("lefse_sort_samples")
        self.on('start', self.start_sort_samples)
        self.on('end', self.end_sort_samples)
        self._memory_increase_step = 40  # 每次重运行增加内存40G by qingchen.zhang @ 20190916

    def start_sort_samples(self):
        self.step.lefse_sort_samples.start()
        self.step.update()

    def end_sort_samples(self):
        self.step.lefse_sort_samples.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        """
        if not self.option("in_otu_table").is_set:
            raise OptionError("输入的丰度文件不能为空", code="32704801")

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["taxa.table.xls", "xls", "各样本物种丰度结果表"],
            ])
        super(LefseSortSampleMgAgent, self).end()

    def set_resource(self):
        """
        设置所需的资源
        """
        self._cpu = 2
        self._memory = "10G"


class LefseSortSampleMgTool(Tool):
    def __init__(self, config):
        super(LefseSortSampleMgTool, self).__init__(config)
        self.logger.info("LefseSortSamplesMg按分组开始处理数据")
        self.abund_table_path = ""

    def cat_samples_percent(self):
        if self.option("group_table").is_set:
            self.abund_table_path = os.path.join(self.output_dir, "taxa.table.xls")
        else:
            raise OptionError("输入的分组文件不能为空", code="32704802")
        sample_groups = {}  # 分组和样本字典
        all_sps = []  # 记录分组文件中的所有样本
        with open(self.option("group_table").path, 'r') as r:
            for l in r:
                if l.startswith('#'):
                    continue
                line = l.strip().split('\t')
                if line[1] not in sample_groups:
                    sample_groups[line[1]] = []
                sample_groups[line[1]].append(line[0])
                all_sps.append(line[0])

        table = pd.read_csv(self.option("in_otu_table").path, sep='\t', index_col=0, chunksize=500000)
        has_header = True
        for one in table:
            one = one[all_sps]
            one.to_csv(self.abund_table_path, sep='\t', mode='a+', header=has_header)
            has_header = False


    def run(self):
        super(LefseSortSampleMgTool, self).run()
        self.cat_samples_percent()
        self.option("out_otu_table").set_path(self.abund_table_path)
        self.end()
