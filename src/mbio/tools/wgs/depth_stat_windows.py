# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.08

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class DepthStatWindowsAgent(Agent):
    """
    软件: depth_stat_windows，处理samtools depth的数据，做覆盖度图
    """
    def __init__(self, parent):
        super(DepthStatWindowsAgent, self).__init__(parent)
        options = [
            {"name": "depth_file", "type": "string"},  # samtools的depth方法的结果文件
            {"name": "step_num", "type": "int", "default": 10000}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("depth_file"):
            raise OptionError("请设置depth_file文件", code="34501901")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(DepthStatWindowsAgent, self).end()


class DepthStatWindowsTool(Tool):
    def __init__(self, config):
        super(DepthStatWindowsTool, self).__init__(config)
        self.depth_stat_windows = "bioinfo/WGS/depth_stat_windows"

    def run_depth_stat_windows(self):
        """
        depth_stat_windows
        """
        sample_name = os.path.basename(self.option("depth_file")).split(".")[0]
        depth_fordraw = sample_name + ".coverage.xls"
        cmd = "{} -i {} -o {} -w {}".format(self.depth_stat_windows, self.option("depth_file"),
                                            self.work_dir + "/" + depth_fordraw, self.option("step_num"))
        command = self.add_command("depth_stat_windows", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("depth_stat_windows运行完成")
        else:
            self.set_error("depth_stat_windows运行失败", code="34501901")
        if os.path.exists(self.output_dir + "/" + depth_fordraw):
            os.remove(self.output_dir + "/" + depth_fordraw)
        os.link(self.work_dir + "/" + depth_fordraw, self.output_dir + "/" + depth_fordraw)

    def run(self):
        super(DepthStatWindowsTool, self).run()
        self.run_depth_stat_windows()
        self.end()
