# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'binbin.zhao'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class VennUpsetAgent(Agent):
    """
    vennupset图
    """
    def __init__(self, parent):
        super(VennUpsetAgent, self).__init__(parent)
        options = [
            {"name": "upset_table", "type": "string"},
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("upset_table"):
            raise OptionError("必须输入upset_table文件")

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
        super(VennUpsetAgent, self).end()


class VennUpsetTool(Tool):
    """
    用于相关性分析
    """
    def __init__(self, config):
        super(VennUpsetTool, self).__init__(config)

    def t_table(self, table_file, new_table):
        """
        转置表格实现
        """
        with open(table_file) as f:
            lines = f.readlines()
            new_table_list = []
            categories_length = len(lines[0].strip().split("\t"))
            for i in range(categories_length):
                new_table_list.append([])
                for line in lines:
                    print line.strip("\n").split("\t")
                    print line.strip("\n").split("\t")[0]
                    try:
                        new_table_list[i].append(line.strip("\n").strip("\r").split("\t")[i])
                    except:
                        new_table_list[i].append("")
            print new_table_list
        with open(new_table, 'w') as w:
            for i in new_table_list:
                w.write("\t".join(i) + "\n")
        return new_table

    def check_group(self):
        with open(self.option("upset_table"), 'r') as f:
            lines = f.readlines()
            categories_length = len(lines[0].strip().split("\t"))
            if categories_length > 10:
                return False
            else:
                return True

    def draw_venn_upset(self):
        """
        开始计算绘图
        """
        self.logger.info("开始绘制venn_upset图")
        if self.check_group():
            self.t_table(self.option("upset_table"), os.path.join(self.output_dir, "venn_upset.txt"))
        else:
            self.logger.info('设置分组过多')
            self.set_error('分组设置最多为十组！', code="32702902")

    def run(self):
        """
        运行
        """
        super(VennUpsetTool, self).run()
        self.draw_venn_upset()  # 对输入数据进行处理，得到new_corr_table.csv
        self.end()
