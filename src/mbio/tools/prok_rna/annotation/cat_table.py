# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class CatTableAgent(Agent):
    """
    CatTable:将fasta文件按行数拆分
    version 1.0
    author: qiuping
    last_modify: 2016.11.15
    """

    def __init__(self, parent):
        super(CatTableAgent, self).__init__(parent)
        options = [
            {"name": "tabledir", "type": "infile", "format": "prok_rna.common_dir"},
            {"name": "tableout", "type": "outfile", "format": "prok_rna.common"},
            {"name": "outname", "type": "string", "default": "merged.xls"},
            {"name": "known", "type": "string", "default": None},
            {"name": "choose_known", "type": "string", "default": None}, #筛选列表
            {"name": "known_format", "type": "string", "default": "{}"},
            {"name": "known_fill", "type": "string", "default": "all"},
            {"name": "head", "type": "bool", "default": False}
        ]
        self.add_option(options)
        self.step.add_steps('cattable')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.cattable.start()
        self.step.update()

    def step_end(self):
        self.step.cattable.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("tabledir").is_set:
            raise OptionError("请传入表格输出结果文件，%s不存在", variables = (self.option("tabledir")), code = "35000301")
        if not self.option("known"):
            pos_num = len(re.findall("{}", self.option("known_format")))
            fill_num = len(self.option("known_fill").split(","))
            if fill_num != fill_num:
                raise OptionError("known_format和known_fill数量不一致", code = "35000302")
        else:
            pass

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '3G'


class CatTableTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(CatTableTool, self).__init__(config)
        self.known_list = []

    def get_list(self):
        '''
        获取需要添加的列名
        '''
        if self.option('choose_known'):
            with open(self.option('choose_known'), 'r') as f:
                self.known_list = [x.strip() for x in f.readlines()]

    def cat_files(self):
        """
        """
        w = open(self.output_dir + "/" + self.option("outname"), 'wb')
        files = os.listdir(self.option('tabledir').prop['path'])
        has_head = self.option("head")
        first_file = True
        for f in files:
            with open(os.path.join(self.option('tabledir').prop['path'], f), 'rb') as r:
                lines = r.readlines()
                if has_head:
                    if first_file:
                        w.writelines(lines)
                    else:
                        w.writelines(lines[1:])
                else:
                    w.writelines(lines)

            first_file = False
            self.option('tableout', self.output_dir +  "/" + self.option("outname"))

        if self.option("known"):
            if self.option("known_format") == "{}" and self.option("known_fill") == "all":
                with open(self.option('known'), 'rb') as r:
                    lines = r.readlines()
                    w.writelines(lines)
            else:
                with open(self.option('known'), 'rb') as r:
                    lines = r.readlines()
                    for line in lines:
                        cols = line.strip().split("\t")
                        if (self.option('choose_known') and cols[0] in self.known_list) or (not self.option('choose_known')):
                            max_col = max([int(i) for i in self.option("known_fill").split(",")])
                            # 表格结构不够填充，则跳过该行
                            if len(cols) < max_col + 1:
                                continue
                            cols_fill = [cols[int(i)] for i in self.option("known_fill").split(",")]
                            w.write(self.option("known_format").format(*tuple(cols_fill)))
        else:
            pass
        w.close()


    def run(self):
        super(CatTableTool, self).run()
        self.get_list()
        self.cat_files()
        self.end()
