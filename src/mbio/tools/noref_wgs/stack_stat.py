# -*- coding: utf-8 -*-
# __author__ = 'Zhaobinbin'
# last_modify: 20181224

import os
import re
import math
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class StackStatAgent(Agent):
    """
    fastq均一化
    """
    def __init__(self, parent=None):
        super(StackStatAgent, self).__init__(parent)
        options = [
            {"name": "snp_vcf", "type": "infile", "format": "dna_gmap.vcf"},
            {"name": "list_file", "type": "infile", "format": "noref_wgs.lists_file"},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("snp_vcf"):
            raise OptionError("snp_vcf is required!", code="35500805")
        if not self.option("list_file"):
            raise OptionError("list_flie is required!", code="35500806")

    def set_resource(self):
        self._cpu = 5   # 因为下面线程为16
        self._memory = "25G"

    def end(self):
        super(StackStatAgent, self).end()


class StackStatTool(Tool):
    def __init__(self, config):
        super(StackStatTool, self).__init__(config)
        self.perl_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/perl"
        self.snp_stat = self.config.PACKAGE_DIR + "/noref_wgs/snp.stat.pl"
        self.variant_qual = self.config.PACKAGE_DIR + "/noref_wgs/variant_qual.pl"
        self.merge = self.config.PACKAGE_DIR + "/noref_wgs/merge.pl"
        self.para_fly = 'program/parafly-r2013-01-21/bin/bin/ParaFly'

    def run_cmd_more(self, cmd_list, cmd_name, cpu):
        """
        将多个cmd命令并行执行
        """
        cmd_file = os.path.join(self.work_dir, "cmd_list_{}.txt".format(cmd_name))
        wrong_cmd = os.path.join(self.work_dir, "failed_cmd_{}.txt".format(cmd_name))
        with open(cmd_file, "w") as f:
            for cmd in cmd_list:
                f.write(cmd + "\n")
        cmd_more = "{} -c {} -CPU {} -failed_cmds {}".format(self.para_fly, cmd_file, cpu, wrong_cmd)
        self.run_cmd(cmd_more, "more_" + cmd_name)


    def run_cmd(self, cmd, cmd_name):
        """
        执行cmd
        """
        self.logger.info(cmd)
        command = self.add_command(cmd_name, cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{}运行完成".format(cmd_name))
        else:
            self.set_error("%s运行失败", code="35500803")

    def run_stack_stat(self):
        cmd_list = []
        cmd_stat = "{} {} -i {} -o {} -m {}".format(self.perl_path, self.snp_stat, self.option("snp_vcf").prop["path"],
                                               self.output_dir + "/" + "snp.stat", self.output_dir + "/" + "snp.matrix")
        cmd_list.append(cmd_stat)
        cmd_qual = "{} {} -i {} -o1 {}".format(self.perl_path, self.variant_qual, self.option("snp_vcf").prop["path"],
                                                    self.output_dir + "/" + "snp.depth")
        cmd_list.append(cmd_qual)
        cmd_merge = "{} {} -i {} -o {}".format(self.perl_path, self.merge, self.option("list_file").prop["path"],
                                               self.output_dir +
                                               "/" + "tag.stat")
        cmd_list.append(cmd_merge)
        self.run_cmd_more(cmd_list, "stack_stat", 3)

    def run(self):
        super(StackStatTool, self).run()
        self.run_stack_stat()
        self.end()
