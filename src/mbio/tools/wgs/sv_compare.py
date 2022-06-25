# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.09

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class SvCompareAgent(Agent):
    """
    SV比较分析
    """
    def __init__(self, parent):
        super(SvCompareAgent, self).__init__(parent)
        options = [
            {"name": "sample1_sv_anno", "type": "infile", 'format': 'bsa.vcf'},  # 样本1的sv注释文件
            {"name": "sample2_sv_anno", "type": "infile", 'format': 'bsa.vcf'},  # 样本2的sv注释文件
            {"name": "is_same", "type": "bool", "default": "true"},  # 选择相同还是不同的结构变异
            {"name": "variation_type", "type": "string"},  # 变异位点类型，"DEL","INV","ITX","CTX","INS"，多个时逗号分隔
            {"name": "variation_len", "type": "string"},  # 变异区域长度,冒号分隔
            {"name": "sample1_support", "type": "string"},  # sample1 reads支持度,冒号分隔
            {"name": "sample2_support", "type": "string"},  # sample2 reads支持度,冒号分隔
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "project_type", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sample1_sv_anno"):
            raise OptionError("请设置sample1的sv注释文件", code="34506901")
        if not self.option("sample2_sv_anno"):
            raise OptionError("请设置sample2的sv注释文件", code="34506902")
        if not self.option("variation_type"):
            raise OptionError("请选择变异位点类型", code="34506903")
        else:
            types = self.option("variation_type").split(",")
            for t in types:
                if t not in ["DEL", "INV", "ITX", "CTX", "INS"]:
                    raise OptionError("变异位点类型%s不在DEL,INV,ITX,CTX,INS内", variables=(t), code="34506904")
        # if self.option("variation_len"):
        #     if not re.search(r".*:.*", self.option("variation_len")):
        #         raise OptionError("变异区域长度{}需用冒号分隔".format(self.option("variation_len")))
        # if self.option("sample1_support"):
        #     if not re.search(r".*:.*", self.option("sample1_support")):
        #         raise OptionError("支持度{}需用冒号分隔".format(self.option("sample1_support")))
        # if self.option("sample2_support"):
        #     if not re.search(r".*:.*", self.option("sample2_support")):
        #         raise OptionError("支持度{}需用冒号分隔".format(self.option("sample2_support")))

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(SvCompareAgent, self).end()


class SvCompareTool(Tool):
    def __init__(self, config):
        super(SvCompareTool, self).__init__(config)
        self.perl_path = "miniconda2/bin/perl"
        self.sv_diff = self.config.PACKAGE_DIR + "/wgs/sv_diff.pl"

    def run_sv_diff(self):
        """
        sv_diff.pl
        """
        if self.option("is_same") == True:
            is_same = "TRUE"
        else:
            is_same = "FALSE"
        cmd = "{} {} -sv1 {} ".format(self.perl_path, self.sv_diff, self.option("sample1_sv_anno").prop['path'])
        cmd += "-sv2 {} -b {} ".format(self.option("sample2_sv_anno").prop['path'], is_same)
        cmd += "-o {} -t {}".format(self.work_dir + "/sv_diff.xls", self.option("variation_type"))
        if self.option("variation_len"):
            cmd += " -l {}".format(self.option("variation_len"))
        if self.option("sample1_support"):
            cmd += " -d1 {}".format(self.option("sample1_support"))
        if self.option("sample2_support"):
            cmd += " -d2 {}".format(self.option("sample2_support"))
        command = self.add_command("sv_diff", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("sv_diff.pl运行完成")
        else:
            self.set_error("sv_diff.pl运行失败", code="34506901")
        os.link(self.work_dir + "/sv_diff.xls", self.output_dir + "/sv_diff.xls")

    def set_db(self):
        """
        导表
        """
        self.logger.info("保存sv比较分析结果到mongo")
        sv_api = self.api.api("wgs.sv")
        if self.option("project_type"):
            sv_api._project_type = self.option("project_type")
        sv_api.add_sg_sv_compare_detail(compare_id=self.option("main_id"), file_path=self.output_dir + "/sv_diff.xls")

    def run(self):
        super(SvCompareTool, self).run()
        self.run_sv_diff()
        self.set_db()
        self.end()
