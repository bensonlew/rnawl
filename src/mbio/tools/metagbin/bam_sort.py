# -*- coding: utf-8 -*-
#__author__ = 'gaohao'@20190729
import os
import re,shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError

class BamSortAgent(Agent):
    """
    将拆分的bam进行sort、index
    """
    def __init__(self, parent):
        super(BamSortAgent, self).__init__(parent)
        options = [
            {"name": "bam", "type": "infile", "format":"align.bwa.bam"}, #输入比对的文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("bam"):
            raise OptionError("必须添加bam文件！", code="")

    def set_resource(self):
        self._cpu = 4
        self._memory = '20G'
        

    def end(self):
        super(BamSortAgent, self).end()


class BamSortTool(Tool):
    def __init__(self, config):
        super(BamSortTool, self).__init__(config)
        self.samtools = "bioinfo/align/samtools-1.4/bin/samtools"
        self.sample_name = os.path.basename(self.option('bam').prop['path']).split('.')[0]

    def run_sort(self):
        """
        对提取的文件进行排序
        :return:
        """
        self.sort_bam = self.output_dir + "/" + self.sample_name + '_sorted.bam'
        for i in os.listdir(self.output_dir):
            os.remove(self.output_dir + "/" + i)
        self.out_bam = self.option('bam').prop['path']
        cmd_sort = '{} sort {} -o {} -@ 4 -m 2G'.format(self.samtools, self.out_bam, self.sort_bam)
        self.logger.info(cmd_sort)
        to_sort = 'to_sort'
        command2 = self.add_command(to_sort, cmd_sort).run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("排序%s运行成功" %to_sort)
        else:
            self.set_error("运行失败")
        self.logger.info("排序运行结束")

    def run_index(self):
        cmd_index = "{} index {}".format(self.samtools, self.sort_bam)
        self.logger.info(cmd_index)
        to_index = 'to_index'
        command2 = self.add_command(to_index, cmd_index).run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("构建索引%s运行成功" % to_index)
        else:
            self.set_error("构建索引运行失败")
        self.logger.info("构建索引运行结束")

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        pass

    def run(self):
        super(BamSortTool, self).run()
        self.run_sort()
        self.run_index()
        self.set_output()
        self.end()