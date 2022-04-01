# -*- coding: utf-8 -*-
#__author__ = 'gaohao'@20190902
import os
import re,shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError

class SampleSortBamAgent(Agent):
    """
    按样品将bam文件合并，再将的bam进行sort、index
    """
    def __init__(self, parent):
        super(SampleSortBamAgent, self).__init__(parent)
        options = [
            {"name": "bam", "type": "string"}, #输入比对的文件
            {"name": "sample", "type": "string"},  # 输入比对的文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("bam"):
            raise OptionError("必须添加bam文件！", code="")

    def set_resource(self):
        self._cpu = 8
        self._memory = '40G'

    def end(self):
        super(SampleSortBamAgent, self).end()


class SampleSortBamTool(Tool):
    def __init__(self, config):
        super(SampleSortBamTool, self).__init__(config)
        self.samtools = "bioinfo/align/samtools-1.4/bin/samtools"
        self.sample_name = self.option('sample')
        self.bam = self.option('bam')

    def run_merge(self):
        """
        按照sample将bam进行merge
        :return:
        """
        if os.path.exists(self.work_dir + '/' + self.sample_name + ".bam"):
            os.remove(self.work_dir + '/' + self.sample_name + ".bam")
        self.merge_bam = self.work_dir + '/' + self.sample_name + ".bam"
        self.logger.info(self.merge_bam)
        cmd_bam = "{} merge --threads 8 {} {}".format(self.samtools, self.merge_bam, self.bam)
        self.logger.info(cmd_bam)
        to_bam = 'merge_bam'
        command1 = self.add_command(to_bam, cmd_bam).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("%s运行完成" %to_bam)
        elif command1.return_code == 1:
            self.add_state('merge.bam文件已经存在', 'memory is low!')
        else:
            self.set_error("%s运行失败", code="")

    def run_sort(self):
        """
        对提取的文件进行排序
        :return:
        """
        self.sort_bam = self.output_dir + "/" + self.sample_name + '.sorted.bam'
        cmd_sort = '{} sort {} -o {} -@ 4 -m 2G'.format(self.samtools, self.merge_bam, self.sort_bam)
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
        super(SampleSortBamTool, self).run()
        self.run_merge()
        self.run_sort()
        self.run_index()
        self.set_output()
        self.end()