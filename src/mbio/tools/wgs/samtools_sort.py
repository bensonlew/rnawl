# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.04

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class SamtoolsSortAgent(Agent):
    """
    软件: samtools
    samtools的sort方法
    """
    def __init__(self, parent):
        super(SamtoolsSortAgent, self).__init__(parent)
        options = [
            {"name": "merged_bam_file", "type": "infile", "format": "align.bwa.bam"},  # bam文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("merged_bam_file").is_set:
            raise OptionError("请设置bam文件", code="34505301")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        super(SamtoolsSortAgent, self).end()


class SamtoolsSortTool(Tool):
    def __init__(self, config):
        super(SamtoolsSortTool, self).__init__(config)
        self.samtools_path = "bioinfo/align/samtools-1.7/samtools"

    def run_samtools_sort(self):
        """
        samtools sort
        """
        for f in os.listdir(self.work_dir):
            if f.endswith(".bam"):
                os.remove(os.path.join(self.work_dir, f))
        sample_name = os.path.basename(self.option("merged_bam_file").prop["path"]).split(".merged.bam")[0]
        bam_name = sample_name + ".sort.bam"
        if os.path.exists(self.work_dir + "/" + bam_name):
            os.remove(self.work_dir + "/" + bam_name)
        cmd = "{} sort -o {} --output-fmt BAM -@ 8 {}".format(self.samtools_path, self.work_dir + "/" + bam_name, self.option("merged_bam_file").prop["path"])
        command = self.add_command("samtools_sort", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("samtools sort完成")
        else:
            i = 0
            while i < 5:
                command = self.add_command("samtools_sort_".format(str(i)), cmd).run()
                self.wait()
                if command.return_code == 0:
                    self.logger.info("samtools sort完成")
                    break
                else:
                    i += 1
                    if i == 5:
                        self.set_error("samtools sort失败", code="34505301")
        if os.path.exists(self.output_dir + "/" + bam_name):
            os.remove(self.output_dir + "/" + bam_name)
        os.link(self.work_dir + "/" + bam_name, self.output_dir + "/" + bam_name)

    def run(self):
        super(SamtoolsSortTool, self).run()
        self.run_samtools_sort()
        self.end()
