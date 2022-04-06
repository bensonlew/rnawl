# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 20190114

import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class SnpStatAgent(Agent):
    """
    ipyrad/stacks 统计snp信息
    """
    def __init__(self, parent=None):
        super(SnpStatAgent, self).__init__(parent)
        options = [
            {"name": "vcf_path", "type": "infile", "format": "noref_wgs.list_file", "required": True},  # ipyrad的data.vcf，stack的populations.snps.vcf
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("vcf_path").is_set:
            raise OptionError("请设置vcf_path", code="35501403")

    def set_resource(self):
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        super(SnpStatAgent, self).end()


class SnpStatTool(Tool):
    def __init__(self, config):
        super(SnpStatTool, self).__init__(config)
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.python_path = "miniconda2/bin/python"
        self.graph_path = self.config.PACKAGE_DIR + "/noref_wgs/snp_depth_density.py"
        self.stat_path = self.config.PACKAGE_DIR + "/noref_wgs/snp.stat.pl"

    def run_snp_stat(self):
        """
        统计snp信息表
        """
        snp_stat = os.path.join(self.output_dir, "snp_stat.xls")
        snp_matrix = os.path.join(self.work_dir, "sample_tag_depth.xls")
        cmd = "{} {} -i {} ".format(self.perl_path, self.stat_path, self.option("vcf_path").prop["path"])
        cmd += "-o {} -m {}".format(snp_stat, snp_matrix)
        command = self.add_command("snp_stat", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("snp信息表统计成功")
        else:
            self.set_error("snp信息表统计失败，请检查", code="35501405")

    def run_snp_depth_density(self):
        """
        统计snp深度和密度分布
        """
        cmd = "{} {} -i {} -o {}".format(self.python_path, self.graph_path, self.option("vcf_path").prop["path"], self.output_dir)
        command = self.add_command("snp_graph", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("统计snp深度和密度分布成功")
        else:
            self.set_error("统计snp深度和密度分布失败，请检查", code="35501406")

    def run(self):
        super(SnpStatTool, self).run()
        self.run_snp_stat()
        self.run_snp_depth_density()
        self.end()
