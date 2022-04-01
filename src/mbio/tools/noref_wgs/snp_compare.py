# -*- coding: utf-8 -*-
# __author__ = 'Zhaobinbin'
# last_modify: 20181218

import os
import re
import math
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class SnpCompareAgent(Agent):
    """
    fastq均一化
    """
    def __init__(self, parent=None):
        super(SnpCompareAgent, self).__init__(parent)
        options = [
            {"name": "vcf_file", "type": "infile", "format": "dna_gmap.vcf"},  # 需要修改infile文件
            {"name": "tag_file", "type": "infile", "format": "noref_wgs.tag_file"},  # 需要修改infile文件
            {"name": "analysis_name", "type": "string"},  # 分析名称
            {"name": "config_file", "type": "infile", "format": "dna_gmap.vcf"},  # infile文件
            {"name": "analysis_model", "type": "string"},

        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("vcf_file"):
            raise OptionError("请输入vcf文件", code="35501111")
        if not self.option("tag_file"):
            raise OptionError("请输入tag文件", code="35501112")
        if not self.option("analysis_name"):
            raise OptionError("请输入分析名称", code="35501113")
        if not self.option("config_file"):
            raise OptionError("请输入config文件", code="35501114")
        if not self.option("analysis_model"):
            raise OptionError("请输入analysis_model文件", code="35501115")

    def set_resource(self):
        self._cpu = 2
        self._memory = "15G"

    def end(self):
        super(SnpCompareAgent, self).end()


class SnpCompareTool(Tool):
    def __init__(self, config):
        super(SnpCompareTool, self).__init__(config)
        self.perl_path = "program/perl/perls/perl-5.24.0/bin/perl"
        self.multiple = self.config.PACKAGE_DIR + "/noref_wgs/noref_SNPselection-Multiple.pl"
        self.single = self.config.PACKAGE_DIR + "/noref_wgs/noref_SNPselection-Single.pl"
        self.stat = self.config.PACKAGE_DIR + "/noref_wgs/compare_stat.py"
        self.Python_path = '/program/Python/bin/python'

    def run_snp_compare(self):

        perl_path = ""
        if self.option("analysis_model") == "single":
            perl_path = self.single
        elif self.option("analysis_model") == "multiple":
            perl_path = self.multiple
        if perl_path == "":
            self.set_error("perl脚本路径不正确", code="35501107")
        cmd = "{} {} -vcf {} -out {} -tag {} -key {} -config {}".format(self.perl_path, perl_path,
                                                                        self.option("vcf_file").prop["path"],
                                                                        self.output_dir,
                                                                        self.option("tag_file").prop["path"],
                                                                        self.option("analysis_name"),
                                                                        self.option("config_file").prop["path"])
        command = self.add_command("snp_compare", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("snp_compare运行完成")
        else:
            self.set_error("snp_compare运行失败", code="35501108")

    def check_vcf_none(self, vcf_path):
        """
        默认vcf不是空的以及是一个样本的分组
        :param vcf_path:
        :return:
        """
        with open(vcf_path, "r") as f:
            num = 0
            for lines in f:
                if not re.match('#', lines):  # 检查vcf不为空
                    num = 1
                    break
            if num == 0:
                return True
            else:
                return False

    def run_snp_compare_stat(self):

        cmd = "{} {} -vcf {} -o {}".format(self.Python_path, self.stat, self.output_dir + "/" +
                                           self.option("analysis_name") + ".vcf", self.output_dir + "/" +
                                           self.option("analysis_name") + ".results")
        command = self.add_command("snp_compare_stat", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("snp_compare_stat运行完成")
        else:
            self.set_error("snp_compare_stat运行失败", code="35501109")

    def run(self):
        super(SnpCompareTool, self).run()
        self.run_snp_compare()
        if not self.check_vcf_none(self.output_dir + "/" + self.option("analysis_name") + ".vcf"):
            self.run_snp_compare_stat()
        self.end()
