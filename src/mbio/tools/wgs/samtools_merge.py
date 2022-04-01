# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.04

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class SamtoolsMergeAgent(Agent):
    """
    软件: samtools
    samtools的merge方法
    """
    def __init__(self, parent):
        super(SamtoolsMergeAgent, self).__init__(parent)
        options = [
            {"name": "bam_list", "type": "string"},  # bam list文件,将bam文件用分号分隔
            {"name": "specimen_id", "type": "string"}  # 输出样本名
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("bam_list"):
            raise OptionError("请设置bam_list", code="34505201")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        super(SamtoolsMergeAgent, self).end()


class SamtoolsMergeTool(Tool):
    def __init__(self, config):
        super(SamtoolsMergeTool, self).__init__(config)
        self.samtools_path = "bioinfo/align/samtools-1.7/samtools"

    def run_samtools_merge(self):
        """
        samtools merge
        """
        self.bam_list = self.option("bam_list").split(";")
        if self.option("specimen_id"):
            sample_name = self.option("specimen_id")
        else:
            sample_name = os.path.basename(self.bam_list[0]).split(".")[0]
        bam_name = sample_name + ".merged.bam"
        cmd = "{} merge -f -c -p -@ 8 --output-fmt BAM {} {}".format(self.samtools_path, self.work_dir + "/" + bam_name, " ".join(self.bam_list))
        command = self.add_command("samtools_merge", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("samtools merge完成")
        else:
            self.set_error("samtools merge失败", code="34505201")
        if os.path.exists(self.output_dir + "/" + bam_name):
            os.remove(self.output_dir + "/" + bam_name)
        os.link(self.work_dir + "/" + bam_name, self.output_dir + "/" + bam_name)

    def run(self):
        super(SamtoolsMergeTool, self).run()
        self.run_samtools_merge()
        self.end()
