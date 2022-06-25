# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.04

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class SamtoolsViewAgent(Agent):
    """
    软件: samtools
    samtools的view方法,将sam文件转为bam文件
    """
    def __init__(self, parent):
        super(SamtoolsViewAgent, self).__init__(parent)
        options = [
            {"name": "sam_file", "type": "infile", "format": "align.bwa.sam"},  # sam文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sam_file").is_set:
            raise OptionError("请设置sam文件", code="34505501")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        super(SamtoolsViewAgent, self).end()


class SamtoolsViewTool(Tool):
    def __init__(self, config):
        super(SamtoolsViewTool, self).__init__(config)
        self.samtools_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/samtools"
        self.samtools_sh_path = "bioinfo/WGS/samtools_view.sh"

    def run_samtools_view(self):
        """
        samtools view
        """
        sample_name = os.path.basename(self.option("sam_file").prop["path"]).split(".sam")[0]
        bam_name = sample_name + ".bam"
        cmd = "{} {} {} {}".format(self.samtools_sh_path, self.samtools_path, self.option("sam_file").prop["path"], self.work_dir + "/" + bam_name)
        command = self.add_command("samtools_view", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("samtools view完成")
        else:
            self.set_error("samtools view失败", code="34505501")
        if os.path.exists(self.output_dir + "/" + bam_name):
            os.remove(self.output_dir + "/" + bam_name)
        os.link(self.work_dir + "/" + bam_name, self.output_dir + "/" + bam_name)

    def run(self):
        super(SamtoolsViewTool, self).run()
        self.run_samtools_view()
        self.end()
