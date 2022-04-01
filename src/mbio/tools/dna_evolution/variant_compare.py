# -*- coding: utf-8 -*-
# __author__ = 'Zhaobinbin'
# modified 2018.05.03

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class VariantCompareAgent(Agent):
    """
    接口“变异位点比较分析Tool”
    """

    def __init__(self, parent):
        super(VariantCompareAgent, self).__init__(parent)
        options = [
            {"name": "variant_compare_config", "type": "infile", "format": "dna_evolution.variant_compare_config"},  # 参数配置文件
            {"name": "filter_recode_vcf", "type": "infile", "format": "sequence.vcf"},
            {"name": "group_table", "type": "infile", "format": "dna_evolution.group_table"}  # 分组信息
        ]
        self.add_option(options)
        self.step.add_steps('variant_compare')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.variant_compare.start()
        self.step.update()

    def step_end(self):
        self.step.variant_compare.finish()
        self.step.update()

    def check_options(self):
        if not self.option("variant_compare_config"):
            raise OptionError("必须输入variant_compare_config参数", code="34902701")
        if not self.option("filter_recode_vcf"):
            raise OptionError("必须输入filter_recode_vcf文件", code="34902702")
        # if not self.option("group_table"):
        #     raise OptionError("必须输入group_table文件")

    def set_resource(self):
        self._cpu = 2
        self._memory = "20G"

    def end(self):
        super(VariantCompareAgent, self).end()


class VariantCompareTool(Tool):

    def __init__(self, config):
        super(VariantCompareTool, self).__init__(config)
        self.perl_path = 'program/perl-5.24.0/bin/perl'
        self.variant_filter_path = self.config.PACKAGE_DIR + "/dna_evolution/variant_filter.pl"
        self.python_path = 'program/Python/bin/python'
        self.chromosome_window_path = self.config.PACKAGE_DIR + "/dna_evolution/step_calc.py"

    def run_variant_compare(self):

        cmd = "{} {} -vcf {} --config {} --out {}".format(self.perl_path, self.variant_filter_path,
                                                          self.option("filter_recode_vcf").prop['path'],
                                                          self.option("variant_compare_config").prop["path"],
                                                          self.output_dir)
        if self.option("group_table").is_set:
            cmd += " --group {}".format(self.option("group_table").prop["path"])
        command = self.add_command("variant_compare", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("variant_compare运行完成")
        else:
            self.set_error("variant_compare运行失败", code="34902701")

        cmd1 = "{} {} -p {} -s {} -o {} -v {}".format(
            self.python_path,
            self.chromosome_window_path,
            os.path.join(self.output_dir, "pop.table"),
            # self.option("pop_table").prop["path"],
            10000,
            self.output_dir,
            "all")
        command1 = self.add_command("chromosome_windows", cmd1).run()
        self.wait()
        if command1.return_code == 0:
            self.logger.info("chromosome_windows运行完成")
        else:
            self.set_error("chromosome_windows运行失败")

    def run(self):
        super(VariantCompareTool, self).run()
        self.run_variant_compare()
        self.end()
