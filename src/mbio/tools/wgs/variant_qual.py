# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180416

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class VariantQualAgent(Agent):
    """
    对snp与indel的质量评估结果，包含质量分布与深度分布
    """
    def __init__(self, parent):
        super(VariantQualAgent, self).__init__(parent)
        options = [
            {"name": "anno_primary_vcf", "type": "string"},
            {"name": "types", "type": "string", "default": 'snp'}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("anno_primary_vcf"):
            raise OptionError("请设置anno_primary_vcf文件", code="34507301")
        if not self.option("types"):
            raise OptionError("请设置types参数--snp/indel", code="34507302")
        else:
            if self.option("types") not in ["snp", "indel"]:
                raise OptionError("%s类型不合法--必须为snp与indel", variables=(self.option("types")), code="34507303")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(VariantQualAgent, self).end()


class VariantQualTool(Tool):
    def __init__(self, config):
        super(VariantQualTool, self).__init__(config)
        self.perl_path = "miniconda2/bin/perl "
        self.script = self.config.PACKAGE_DIR + "/wgs/variant_qual.pl"

    def run_tool(self):
        """
        perl variant_qual.pl -i indel.anno.primary.vcf -o1 indel.depth -o2 indel.GQ

        perl variant_qual.pl -i snp.anno.primary.vcf -o1 snp.depth -o2 snp.GQ
        """
        cmd = "{}{}".format(self.perl_path, self.script)
        if self.option("types") == 'snp':
            cmd += " -i {} -o1 {} -o2 {}"\
                .format(self.option("anno_primary_vcf"), os.path.join(self.output_dir, "snp.depth"),
                        os.path.join(self.output_dir, "snp.GQ"))
        else:
            cmd += " -i {} -o1 {} -o2 {}" \
                .format(self.option("anno_primary_vcf"), os.path.join(self.output_dir, "indel.depth"),
                        os.path.join(self.output_dir, "indel.GQ"))
        command = self.add_command("script", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("脚本运行完成!")
        else:
            self.set_error("脚本运行失败", code="34507301")

    def run(self):
        super(VariantQualTool, self).run()
        self.run_tool()
        self.end()
