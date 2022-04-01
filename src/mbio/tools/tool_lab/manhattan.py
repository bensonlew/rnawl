# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'binbin.zhao'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class ManhattanAgent(Agent):
    """
    CorHeatmapAgent:用于生成之间的correlation
    """
    def __init__(self, parent):
        super(ManhattanAgent, self).__init__(parent)
        options = [
            {"name": "snp_table", "type": "infile", "format": "tool_lab.snp_table"},
            {"name": "guides2", "type": "string", "default": "False"},  # 第二条阈值线是否画, 如果话传yes。
            {"name": "p_threshold1", "type": "float", "default": 1e-05},  # 当不画的时候直接传0过来即可。
            {"name": "p_threshold2", "type": "float", "default": 5e-08},
            {"name": "color", "type": "string"},
            {"name": "re_axis", "type": "string", "default": -1},
            {"name": "one_chr", "type": "int", "default": -1},
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("snp_table"):
            raise OptionError("必须输入snp_table文件")
        if not self.option("color"):
            raise OptionError("必须输入散点颜色")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        """
        计算结束
        """
        super(ManhattanAgent, self).end()


class ManhattanTool(Tool):
    """
    用于相关性分析
    """
    def __init__(self, config):
        super(ManhattanTool, self).__init__(config)
        self.r_path = '/program/R-3.3.1/bin/Rscript'
        self.script_path = self.config.PACKAGE_DIR + "/tool_lab/qqman.R"


    def draw_manhattan(self):
        """
        开始计算绘图
        """
        self.logger.info(self.option("color"))
        cmd = "{} {} -s {}  -o {} -c {} -l {} -y {} -g {}".format(
            self.r_path,
            self.script_path,
            self.option("snp_table").prop['path'],
            self.output_dir,
            self.option("color"),
            self.option("p_threshold1"),
            self.option("guides2"),
            self.option("p_threshold2")
            )
        self.logger.info("manhattan绘制完成")
        if self.option("re_axis") != "-1":
        # if self.option("re_axis").encode("utf-8") != "-1":
            cmd += " -b {}".format(self.option("re_axis"))
        if self.option("one_chr") != -1:
            cmd += " -p {}".format(self.option("one_chr"))
        self.logger.info(cmd)
        command = self.add_command("manhattan", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("manhattan绘制完成")
        else:
            self.set_error("manhattan绘制失败")


    def run(self):
        """
        运行
        """
        super(ManhattanTool, self).run()
        self.draw_manhattan() # 对输入数据进行处理，得到new_corr_table.csv
        self.end()