# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190320

import os
import re
import math
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class SsrMisaAgent(Agent):
    """
    将每个样本的ssr结果进行合并，进行misa和统计，得到每个样本的ssr.result和分类统计结果
    """
    def __init__(self, parent):
        super(SsrMisaAgent, self).__init__(parent)
        options = [
            {"name": "ssr_list", "type": "infile", "format": "wgs_v2.bam_list", "required": True},  # 样本对应的ssr.result结果文件list
            {"name": "sample_name", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("ssr_list").is_set:
            raise OptionError("请设置样本的ssr.result结果文件list")

    def set_resource(self):
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(SsrMisaAgent, self).end()


class SsrMisaTool(Tool):
    def __init__(self, config):
        super(SsrMisaTool, self).__init__(config)
        self.perl_path = "miniconda2/bin/perl"
        self.python_path = "miniconda2/bin/python"
        self.ssr_misa_path = self.config.PACKAGE_DIR + "/wgs_v2/sample_ssr_misa.pl"
        self.ssr_stat_path = self.config.PACKAGE_DIR + "/wgs_v2/ssr_misa_stat.py"

    def merge_all_ssr_result(self):
        """
        合并所有的ssr.result,对低值的进行过滤
        """
        self.ssr_result = os.path.join(self.output_dir, self.sample_name + ".ssr.result")
        with open(self.option("ssr_list").prop["path"], "r") as f, open(self.ssr_result, "w") as w:
            w.write("#Chr\tPos\tSSRbase\tUnit\tDep\tTotaldep\n")
            for line in f:
                item = line.strip().split("\t")
                with open(item[1], "r") as f1:
                    lines = f1.readlines()
                    for line1 in lines[1:]:
                        item1 = line1.strip().split("\t")
                        unit = item1[3].split(":")
                        dep = item1[4].split(":")
                        unit_, dep_ = [], []
                        if len(unit) > 1:
                            for i in range(len(unit)):
                                if int(dep[i]) < 2:
                                    continue
                                if len(item1[2]) * int(unit[i]) < 8:
                                    continue
                                unit_.append(unit[i])
                                dep_.append(dep[i])
                            if unit_:
                                new = [item1[0], item1[1], item1[2], ":".join(unit_), ":".join(dep_), item1[5]]
                                w.write("\t".join(new) + "\n")
                        else:
                            if len(item1[2]) * int(unit[0]) >= 8:
                                w.write(line1)

    def run_sample_ssr_misa(self):
        """
        sample_ssr_misa.pl:对合并的ssr.result进行misa，得到ssr.misa
        """
        self.ssr_misa = os.path.join(self.output_dir, self.sample_name + ".ssr.misa")
        self.logger.info(self.ssr_misa_path)
        cmd = "{} {} -i {} -o {}".format(self.perl_path, self.ssr_misa_path, self.ssr_result, self.ssr_misa)
        self.logger.info(cmd)
        command = self.add_command("sample_ssr_misa", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("对样本进行misa成功")
        else:
            self.set_error("对样本进行misa失败，请检查")

    def run_sample_ssr_stat(self):
        """
        对misa的结果ssr.misa进行统计，得到SSR类型分布进行统计
        """
        ssr_stat = os.path.join(self.output_dir, self.sample_name + ".ssr.stat.xls")
        cmd = "{} {} -i {} -o {} -s {}".format(self.python_path, self.ssr_stat_path, self.ssr_misa, ssr_stat, self.option("sample_name"))
        command = self.add_command("ssr_stat", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("SSR统计成功")
        else:
            self.set_error("SSR统计失败，请检查")

    def run(self):
        super(SsrMisaTool, self).run()
        self.sample_name = self.option("sample_name") if self.option("sample_name") else "merge"
        self.merge_all_ssr_result()
        self.run_sample_ssr_misa()
        self.run_sample_ssr_stat()
        self.end()
