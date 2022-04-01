# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'binbin.zhao'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class CircosAgent(Agent):
    """
    vennupset图
    """
    def __init__(self, parent):
        super(CircosAgent, self).__init__(parent)
        options = [
            {"name": "circos_file", "type": "infile", "format": "tool_lab.upset_table"},
            {"name": "circos_type", "type": "string"},
            {"name": "index", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("circos_file"):
            raise OptionError("请设置circos_file")
        if not self.option("circos_type"):
            raise OptionError("请设置circos_type")

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
        super(CircosAgent, self).end()


class CircosTool(Tool):
    """
    用于相关性分析
    """
    def __init__(self, config):
        super(CircosTool, self).__init__(config)

    def draw_circos_chart(self):
        """
        开始计算绘图
        """
        with open(self.option("circos_file").prop['path']) as f, open(os.path.join(self.output_dir, "circos_chart.json"
                                                              + self.option("index")), "w") as w:
            lines = f.readlines()
            for line in lines:
                item = line.strip().split("\t")
                self.logger.info(item)
                w.write("{chr:" + item[0].strip() + ",start:" + item[1].strip() + ",end:" + item[2].strip() + ",value:"
                        + item[3].strip() + ",group:" + item[4].strip() + "},\n")

    def draw_circos_highlight(self):
        """
        转化highlight的数据结构
        """
        with open(self.option("circos_file").prop['path']) as f, open(os.path.join(self.output_dir, "circos_highlight.json"
                                                              + self.option("index")), "w") as w:
            lines = f.readlines()
            for line in lines:
                item = line.strip().split("\t")
                w.write("{chr:" + item[0].strip() + ",start:" + item[1].strip() + ",end:" + item[2].strip() + ",color:"
                        + item[3].strip() + "},\n")


    def draw_circos_chr(self):
        """
        转化highlight的数据结构
        """
        with open(self.option("circos_file").prop['path']) as f, open(os.path.join(self.output_dir, "circos_chr.json"
                                                              + self.option("index")), "w") as w:
            lines = f.readlines()
            for line in lines:
                item = line.strip().split("\t")
                w.write("[" + item[0].strip() + "," + item[1].strip() + "]," + "\n")

    def run(self):
        """
        运行
        """
        super(CircosTool, self).run()
        if self.option("circos_type") == "c_highlight":
            self.draw_circos_highlight()
        elif self.option("circos_type") == "c_chr":
            self.draw_circos_chr()
        else:
            self.draw_circos_chart()  # 对输入数据进行处理，得到new_corr_table.csv
        self.end()
