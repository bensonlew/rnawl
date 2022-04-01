# -*- coding: utf-8 -*-
# __author__ = 'wentian.liu'
# modified 2019.02.20

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class BcftoolsConvertAgent(Agent):
    """
    工具：bcftools convert
    """
    def __init__(self, parent):
        super(BcftoolsConvertAgent, self).__init__(parent)
        options = [
            {"name": "bcf_file", "type": "infile", "format": "wgs_v2.bcf"},
            {"name": "vcf_file", "type": "outfile", "format": "sequence.vcf"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("bcf_file").is_set:
            raise OptionError("请设置bcf_file")

    def set_resource(self):
        self._cpu = 2
        self._memory = "20G"

    def end(self):
        super(BcftoolsConvertAgent, self).end()


class BcftoolsConvertTool(Tool):
    def __init__(self, config):
        super(BcftoolsConvertTool, self).__init__(config)
        self.bcftools = "bioinfo/seq/bcftools-1.7/bcftools"

    def run_bcftools_convert(self):
        """
        bcftools convert -O v -o $dOut/pop.sort.sv.vcf $dOut/pop.nosort.sv.bcf
        """
        output_path = os.path.join(self.output_dir, "pop.sort.sv.vcf")
        cmd = "{} convert -O v -o {} {}".format(self.bcftools, output_path, self.option("bcf_file").prop["path"])
        command = self.add_command("bcftools_convert", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("bcftools_convert运行成功")
        else:
            self.set_error("bcftools_convert运行失败")
        self.option("vcf_file", output_path)

    def run(self):
        super(BcftoolsConvertTool, self).run()
        self.run_bcftools_convert()
        self.end()
