# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180416

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class VcfStatAgent(Agent):
    """
    对snp与indel注释后的vcf进行分析，snp.anno.primary.vcf/indel.anno.primary.vcf得到功效统计与功能信息统计
    """
    def __init__(self, parent):
        super(VcfStatAgent, self).__init__(parent)
        options = [
            {"name": "anno_primary_vcf", "type": "string"},
            {"name": "types", "type": "string", "default": 'snp'}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("anno_primary_vcf"):
            raise OptionError("请设置anno_primary_vcf文件", code="34507501")
        if not self.option("types"):
            raise OptionError("请设置types参数--snp/indel", code="34507502")
        else:
            if self.option("types") not in ["snp", "indel"]:
                raise OptionError("%s类型不合法--必须为snp与indel", variables=(self.option("types")), code="34507503")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(VcfStatAgent, self).end()


class VcfStatTool(Tool):
    def __init__(self, config):
        super(VcfStatTool, self).__init__(config)
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl "
        self.script = self.config.PACKAGE_DIR + "/wgs/vcfstat.pl"

    def run_tool(self):
        """
        perl /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/wgs/vcfstat.pl
         -i /mnt/ilustre/users/sanger-dev/workspace/20180416/Single_snp_eff_20180416/Snpeff/output/snp.anno.primary.vcf
         -o test
        """
        cmd = "{}{} -i {} -o {}".format(self.perl_path, self.script, self.option("anno_primary_vcf"), self.work_dir)
        command = self.add_command("script", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("脚本运行完成!")
        else:
            self.set_error("脚本运行失败", code="34507501")

        self.logger.info("开始设置结果目录！")
        if self.option("types") == 'snp':
            if os.path.exists(self.output_dir + "/snp.stat"):
                os.remove(self.output_dir + "/snp.stat")
            os.link(self.work_dir + "/snp.stat", self.output_dir + "/snp.stat")
        else:
            if os.path.exists(self.output_dir + "/indel.stat"):
                os.remove(self.output_dir + "/indel.stat")
            os.link(self.work_dir + "/indel.stat", self.output_dir + "/indel.stat")
        self.logger.info("设置结果目录成功！")

    def run(self):
        super(VcfStatTool, self).run()
        self.run_tool()
        self.end()
