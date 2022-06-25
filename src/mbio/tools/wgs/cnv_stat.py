# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.10

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class CnvStatAgent(Agent):
    """
    cnv 变异类型统计,变异长度统计
    """
    def __init__(self, parent):
        super(CnvStatAgent, self).__init__(parent)
        options = [
            {"name": "cnv_anno", "type": "string"},  # cnv注释文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("cnv_anno"):
            raise OptionError("请设置cnv注释文件", code="34501601")
        if not os.path.exists(self.option("cnv_anno")):
            raise OptionError("cnv注释文件:%s不存在，请检查",variables=(self.option("cnv_anno")), code="34501602")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(CnvStatAgent, self).end()


class CnvStatTool(Tool):
    def __init__(self, config):
        super(CnvStatTool, self).__init__(config)
        self.perl_path = "miniconda2/bin/perl"
        self.cnv_stat = self.config.PACKAGE_DIR + "/wgs/cnv_stat.pl"

    def run_cnv_stat(self):
        """
        cnv_stat.pl cnv变异类型和长度统计
        """
        self.sample_name = os.path.basename(self.option("cnv_anno")).split(".cnv.anno")[0]
        cmd = "{} {} -i {} -o {}".format(self.perl_path, self.cnv_stat, self.option("cnv_anno"), self.work_dir + "/" + self.sample_name)
        command = self.add_command("cnv_stat", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("cnv_stat.pl运行完成")
        else:
            self.set_error("cnv_stat.pl运行失败", code="34501601")

    def set_output(self):
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        os.link(self.work_dir + "/" + self.sample_name + ".cnv.stat", self.output_dir + "/" + self.sample_name + ".cnv.stat.xls")
        os.link(self.work_dir + "/" + self.sample_name + ".cnv.length", self.output_dir + "/" + self.sample_name + ".cnv.length.xls")

    def run(self):
        super(CnvStatTool, self).run()
        self.run_cnv_stat()
        self.set_output()
        self.end()
