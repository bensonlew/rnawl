# -*- coding: utf-8 -*-
# __author__ = 'liuwentian'
# modified 2018.06.28

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class GetGroupingResultAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(GetGroupingResultAgent, self).__init__(parent)
        options = [
            {"name": "pop_marker_detail", "type": "infile", "format": "dna_gmap.marker"},  # pop.filterd.detail.info文件
            {"name": "total_lg", "type": "infile", "format": "dna_gmap.lg"},  # total.lg文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("total_lg").is_set:
            raise OptionError("请设置total.lg文件", code="34800201")
        if not self.option("pop_marker_detail").is_set:
            raise OptionError("请设置pop.filterd.detail.info文件", code="34800202")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        super(GetGroupingResultAgent, self).end()


class GetGroupingResultTool(Tool):
    def __init__(self, config):
        super(GetGroupingResultTool, self).__init__(config)
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl'
        self.lgmarker_bin = self.config.PACKAGE_DIR + "/dna_gmap/lgmarker_stat.pl"

    def run_marker_stat_bin(self):
        """
        进行标记信息统计
        """
        cmd = "{} {} -lg {}".format(self.perl_path, self.lgmarker_bin, self.option("total_lg").prop["path"])
        cmd += " -info {} -out {}".format(self.option("pop_marker_detail").prop["path"], self.output_dir + "/total.marker.stat.xls")
        command = self.add_command("lgmarker_stat", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("lgmarker_stat完成")
        else:
            self.set_error("lgmarker_stat失败", code="34800201")

    def run(self):
        super(GetGroupingResultTool, self).run()
        self.run_marker_stat_bin()
        self.end()
