# -*- coding: utf-8 -*-
# __author__ = 'Zhaobinbin'
# modified 2019.03.19

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class VariantCompareAgent(Agent):
    """
    接口“变异位点比较分析Tool, 导表在tool中完成”
    """

    def __init__(self, parent):
        super(VariantCompareAgent, self).__init__(parent)
        options = [
            {"name": "variant_compare_config", "type": "infile",
                "format": "wgs_v2.variant_compare_config"},  # 参数配置文件
            {"name": "filter_recode_vcf", "type": "infile", "format": "wgs_v2.vcf"},
            {"name": "group_table", "type": "infile",
                "format": "wgs_v2.group_table"},  # 分组信息
            {"name": "name", "type": "string"},  # 传入文件的名称
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
            raise OptionError("必须输入variant_compare_config参数")
        if not self.option("filter_recode_vcf"):
            raise OptionError("必须输入filter_recode_vcf文件")
        if not self.option("group_table"):
            raise OptionError("必须输入group_table文件")
        if not self.option("name"):
            raise OptionError("必须输入生成文件名称")

    def set_resource(self):
        self._cpu = 2
        self._memory = "20G"

    def end(self):
        super(VariantCompareAgent, self).end()


class VariantCompareTool(Tool):

    def __init__(self, config):
        super(VariantCompareTool, self).__init__(config)
        self.python = '/miniconda2/bin/python'
        self.variant_filter_path = self.config.PACKAGE_DIR + \
            "/wgs_v2/variant_compare_filter.py"

    def run_variant_compare(self):
        self.checkconfig()
        cmd = "{} {} -i {} -c {} -o {} -g {} -n {}".format(
            self.python,
            self.variant_filter_path,
            self.option("filter_recode_vcf").prop['path'],
            self.option("variant_compare_config").prop["path"],
            self.output_dir,
            self.option("group_table").prop["path"], self.option("name"))
        command = self.add_command("variant_compare", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("variant_compare运行完成")
        else:
            self.set_error("variant_compare运行失败")

    def checkconfig(self):
        with open(self.option("variant_compare_config").prop["path"], 'r') as f:
            lines = f.readlines()[3:]
            for line in lines:
                if re.match("Region", line):
                    genome_region = line.strip().split("=")[1]
                    if genome_region != "all":
                        chr_region = genome_region.strip().split(":")[0].lower()
                        if not re.match("^chr.*", chr_region) and not re.match("^sca.*", chr_region):
                            self.set_error("基因组区域：{}不合法，请确保染色体为chr或者sca".format(chr_region))

    def run(self):
        super(VariantCompareTool, self).run()
        self.run_variant_compare()
        self.end()
