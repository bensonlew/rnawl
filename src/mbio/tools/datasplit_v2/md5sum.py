# -*- coding: utf-8 -*-
# __author__ = 'zengjing'

"""md5校验 """
import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class Md5sumAgent(Agent):
    """
    md5校验
    """
    def __init__(self, parent=None):
        super(Md5sumAgent, self).__init__(parent)
        options = [
            {'name': 'fastq_dir', 'type': "string"},  # 进行校验的文件夹
        ]
        self.add_option(options)
        # self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_options(self):
        if not self.option('fastq_dir'):
            raise OptionError('必须输入fastq_dir')
        return True

    def set_resource(self):
        """
        设置所需要的资源
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(Md5sumAgent, self).end()


class Md5sumTool(Tool):
    def __init__(self, config):
        super(Md5sumTool, self).__init__(config)

    def run_md5sum(self):
        """
        生成md5校验码
        """
        fastq_dir = os.path.join(self.option("fastq_dir"))
        self.fq_list = []
        md5_info = {}
        for f in os.listdir(fastq_dir):
            fq_path = os.path.join(fastq_dir, f)
            if os.path.isdir(fq_path):
                continue
            md5 = os.popen("md5sum {}".format(fq_path)).readlines()[0].split(" ")[0]
            md5_info[f] = md5
        with open(os.path.join(fastq_dir, "md5sum.txt"), "wb") as w:
            for f in md5_info.keys():
                w.write(md5_info[f] + "  " + f + "\n")

    def run(self):
        super(Md5sumTool, self).run()
        self.run_md5sum()
        self.end()
