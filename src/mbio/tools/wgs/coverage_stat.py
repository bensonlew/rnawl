# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.08

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class CoverageStatAgent(Agent):
    """
    map_stat.pl
    统计bam文件的depth、insert、stat
    """
    def __init__(self, parent):
        super(CoverageStatAgent, self).__init__(parent)
        options = [
            {"name": "map_stat", "type": "string"},  # samtools的stats结果文件
            {"name": "mkdup_metric", "type": "string"},  # picard的mkdup的结果metric文件
            {"name": "ref_dict", "type": "string"},  # ref.ref_dict
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("map_stat"):
            raise OptionError("请设置bam文件", code="34501801")
        if not self.option("ref_dict"):
            raise OptionError("请设置ref.dict文件", code="34501802")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(CoverageStatAgent, self).end()


class CoverageStatTool(Tool):
    def __init__(self, config):
        super(CoverageStatTool, self).__init__(config)
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.map_stat = self.config.PACKAGE_DIR + "/wgs/map_stat.pl"

    def run_map_stat(self):
        """
        map_stat
        """
        self.sample_name = os.path.basename(self.option("map_stat")).split(".")[0]
        cmd = "{} {} -b {} ".format(self.perl_path, self.map_stat, self.option("map_stat"))
        cmd += "-i {} ".format(os.path.join(self.work_dir, self.sample_name + ".insert"))
        cmd += "-c {} ".format(os.path.join(self.work_dir, self.sample_name + ".depth"))
        cmd += "-d {} -o {} -k {}".format(self.option("ref_dict"), os.path.join(self.work_dir, self.sample_name + ".result.stat"), self.sample_name)
        if self.option("mkdup_metric"):
            cmd += " -m {}".format(self.option("mkdup_metric"))
        command = self.add_command("map_stat", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("map_stat.pl运行完成")
        else:
            self.set_error("map_stat.pl运行失败", code="34501801")

    def set_output(self):
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        os.link(os.path.join(self.work_dir, self.sample_name + ".insert"), os.path.join(self.output_dir, self.sample_name + ".insert.xls"))
        os.link(os.path.join(self.work_dir, self.sample_name + ".depth"), os.path.join(self.output_dir, self.sample_name + ".depth.xls"))
        os.link(os.path.join(self.work_dir, self.sample_name + ".result.stat"), os.path.join(self.output_dir, self.sample_name + ".result.stat.xls"))

    def run(self):
        super(CoverageStatTool, self).run()
        self.run_map_stat()
        self.set_output()
        self.end()
