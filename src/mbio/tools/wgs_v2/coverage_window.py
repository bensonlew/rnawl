# -*- coding: utf-8 -*-
# __author__ = 'binbin zhao'
# modified 2019.03.25

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import datetime

starttime = datetime.datetime.now()
class CoverageWindowAgent(Agent):
    """
    软件: depth_stat_window，处理samtools depth的数据，做覆盖度图
    """
    def __init__(self, parent):
        super(CoverageWindowAgent, self).__init__(parent)
        options = [
            {"name": "bam_file", "type": "infile", "format": "align.bwa.bam"},
            {"name": "step_num", "type": "int", "default": 200000}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("bam_file"):
            raise OptionError("请设置depth_file路径", code="34902601")
        if not self.option("step_num") in [5000, 10000, 50000, 100000, 200000, 500000]:
            raise OptionError("step_num 只能是5000，10000，50000,100000，200000，500000", code="34902602")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(CoverageWindowAgent, self).end()


class CoverageWindowTool(Tool):
    def __init__(self, config):
        super(CoverageWindowTool, self).__init__(config)
        self.depth_stat_window = "bioinfo/WGS/depth_stat_windows"
        self.samtools_sh_path = "bioinfo/WGS/samtools_depth.sh"
        self.samtools_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/samtools"

    def run_depth_stat_window(self):
        """
        depth_stat_window
        """
        sample_name = os.path.basename(self.option("bam_file").prop["path"]).split(".")[0]
        depth_fordraw = sample_name + ".coverage.xls"
        cmd = "{} -i {} -o {} -w {}".format(self.depth_stat_window, self.work_dir + "/" + sample_name + ".coverage",
                                            self.work_dir + "/" + depth_fordraw, self.option("step_num"))
        command = self.add_command("depth_stat_window", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("depth_stat_window运行完成")
        else:
            self.set_error("depth_stat_window运行失败")
        if os.path.exists(self.output_dir + "/" + depth_fordraw):
            os.remove(self.output_dir + "/" + depth_fordraw)
        os.link(self.work_dir + "/" + depth_fordraw, self.output_dir + "/" + depth_fordraw)

    def run_samtools_depth(self):
        """
        先把bam文件转化为run_depth_stat_window 传入文件
        :return:
        """
        sample_name = os.path.basename(self.option("bam_file").prop["path"]).split(".")[0]
        depth_file = sample_name + ".coverage"
        cmd = "{} {} {} {}".format(self.samtools_sh_path, self.samtools_path, self.option("bam_file").prop["path"],
                                   self.work_dir + "/" + depth_file)
        command = self.add_command("samtools_depth", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("samtools depth完成")
        else:
            self.set_error("samtools depth失败")

    def run(self):
        super(CoverageWindowTool, self).run()
        self.run_samtools_depth()
        self.run_depth_stat_window()
        self.end()
