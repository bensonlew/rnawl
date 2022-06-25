# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.05.16

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class SsrCompareStatAgent(Agent):
    """
    SSR比较分析统计
    """
    def __init__(self, parent):
        super(SsrCompareStatAgent, self).__init__(parent)
        options = [
            {"name": "sample1_result", "type": "string"},  # 样本1的result文件
            {"name": "sample2_result", "type": "string"},  # 样本2的result文件
            {"name": "is_same", "type": "bool"},  # 样本1和样本2之间相同或不同的SSR位点
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sample1_result"):
            raise OptionError("请设置样本1的result文件", code="34506201")
        if not self.option("sample2_result"):
            raise OptionError("请设置样本2的result文件", code="34506202")
        if not os.path.exists(self.option("sample1_result")):
            raise OptionError("文件:%s不存在，请检查", ariables=(self.option("sample1_result")), code="34506203")
        if not os.path.exists(self.option("sample2_result")):
            raise OptionError("文件:%s不存在，请检查", variables=(self.option("sample2_result")), code="34506204")

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(SsrCompareStatAgent, self).end()


class SsrCompareStatTool(Tool):
    def __init__(self, config):
        super(SsrCompareStatTool, self).__init__(config)
        self.perl_path = "miniconda2/bin/perl"
        self.ssr_stat = self.config.PACKAGE_DIR + "/wgs/ssr-stat.pl"
        self.ssr_table = self.config.PACKAGE_DIR + "/wgs/ssr2table.pl"

    def run_ssr_stat(self):
        """
        ssr-stat.pl
        """
        type = "same" if self.option("is_same") else "diff"
        cmd = "{} {} -type {} -i {} -k {} -o {}".format(self.perl_path, self.ssr_stat, type, self.option("sample1_result"),
                                                        self.option("sample2_result"), self.work_dir)
        command = self.add_command("ssr_stat", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("ssr-stat.pl运行完成")
        else:
            self.set_error("ssr-stat.pl运行失败", code="34506201")
        ssr_detail = "ssr." + type + ".result.xls"
        if os.path.exists(os.path.join(self.output_dir, "ssr_detail.xls")):
            os.remove(os.path.join(self.output_dir, "ssr_detail.xls"))
        os.link(os.path.join(self.work_dir, ssr_detail), os.path.join(self.output_dir, "ssr_detail.xls"))

    def run_ssr_table(self):
        """
        ssr2table.pl
        """
        stat_file = os.path.join(self.output_dir, "ssr_detail.xls")
        cmd = "{} {} -i {} -o {}".format(self.perl_path, self.ssr_table, stat_file, os.path.join(self.output_dir, "ssr_stat.xls"))
        command = self.add_command("ssr_table", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("ssr2table.pl运行完成")
        else:
            self.set_error("ssr2table.pl运行失败", code="34506202")

    def run(self):
        super(SsrCompareStatTool, self).run()
        self.run_ssr_stat()
        self.run_ssr_table()
        self.end()
