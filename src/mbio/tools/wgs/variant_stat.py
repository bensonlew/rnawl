# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180416

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class VariantStatAgent(Agent):
    """
    对snp与indel注释后的vcf进行分析分析统计snp与indel位点信息统计
    """
    def __init__(self, parent):
        super(VariantStatAgent, self).__init__(parent)
        options = [
            {"name": "anno_primary_vcf", "type": "string"},
            {"name": "types", "type": "string", "default": 'snp'},
            {"name": "need_mutation_distribution", "type": "bool", "default": False}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("anno_primary_vcf"):
            raise OptionError("请设置anno_primary_vcf文件", code="34507401")
        if not self.option("types"):
            raise OptionError("请设置types参数--snp/indel", code="34507402")
        else:
            if self.option("types") not in ["snp", "indel"]:
                raise OptionError("%s类型不合法--必须为snp与indel", variables=(self.option("types")), code="34507403")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(VariantStatAgent, self).end()


class VariantStatTool(Tool):
    def __init__(self, config):
        super(VariantStatTool, self).__init__(config)
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl "
        self.python_path = "miniconda2/bin/python "
        self.script1 = self.config.PACKAGE_DIR + "/wgs/snp.stat.pl"
        self.script2 = self.config.PACKAGE_DIR + "/wgs/indel.stat.pl"
        self.script3 = self.config.PACKAGE_DIR + "/wgs_v2/snp_mutation_distribution.py"

    def run_tool(self):
        """
        perl snp.stat.pl -i snp.anno.primary.vcf -o snp.stat -m snp.matrix
        perl indel.stat.pl -i indel.anno.primary.vcf -o indel.stat -m indel.matrix -l indel.len
        """
        cmd = "{}".format(self.perl_path)
        if self.option("types") == 'snp':
            cmd += "{} -o {} -m {}".format(self.script1, os.path.join(self.output_dir, "snp.stat"),
                                           os.path.join(self.output_dir, "snp.matrix"))
        else:
            cmd += "{} -o {} -m {} -l {}"\
                .format(self.script2, os.path.join(self.output_dir, "indel.stat"),
                        os.path.join(self.output_dir, "indel.matrix"), os.path.join(self.output_dir, "indel.len"))
        cmd += " -i {}".format(self.option("anno_primary_vcf"))
        command = self.add_command("script", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("脚本运行完成!")
        else:
            self.set_error("脚本运行失败", code="34507401")

    def run_mutation_distribution(self):
        """
        python snp_mutation_distribution.py -i -t -o
        indel_gene_distribution.txt
        snp_type_distribution.txt
        :return:
        """
        cmd = "{} {} -i {} -t {} -o {}".format(self.python_path, self.script3, self.option("anno_primary_vcf"),
                                               self.option("types"), self.output_dir)
        command = self.add_command("mutation_distribution", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("脚本运行完成!")
        else:
            self.set_error("脚本运行失败")

    def run(self):
        super(VariantStatTool, self).run()
        self.run_tool()
        if self.option("need_mutation_distribution"):
            self.run_mutation_distribution()
        self.end()
