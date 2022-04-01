# !usr/bin/python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import os


class BubbleAgent(Agent):
    """
    气泡图
    """
    def __init__(self, parent):
        super(BubbleAgent, self).__init__(parent)
        options = [
            {"name": "tooltable", "type": "infile", "format": "tool_lab.table"},
            {"name": "x_axis", "type": "string"},
            {"name": "y_axis", "type": "string"},
            {"name": "size", "type": "string"},
            {"name": "group", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        """
        if not self.option("tooltable"):
            raise OptionError('必须设置作图数据')
        if not self.option("x_axis"):
            raise OptionError('必须设置X轴数据')
        if not self.option("x_axis"):
            raise OptionError('必须设置Y轴数据')
        if not self.option("size"):
            raise OptionError('必须设置气泡大小')

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "筛选气泡图数据"],
            ["./bubble.csv", "csv", "气泡图数据"]
        ])
        super(BubbleAgent, self).end()


class BubbleTool(Tool):
    """
    气泡图
    """
    def __init__(self, config):
        super(BubbleTool, self).__init__(config)
        self._version = 'v2.1-20140214'
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')

    def run(self):
        """
        运行
        """
        super(BubbleTool, self).run()
        self.get_table()
        self.set_output()
        # self.set_db()
        self.end()

    def get_table(self):
        """
        提取用于画图的数据
        """
        table = self.option("tooltable").prop['path']
        df = pd.read_table(table, header=0)
        x = self.option("x_axis")
        y = self.option("y_axis")
        s = self.option("size")
        if self.option("group"):
            g = self.option("group")
            new_df = df[[x, y, s, g]]
            new_df.columns = ['x', 'y', 'size', 'group']
        else:
            new_df = df[[x, y, s]]
            new_df.columns = ['x', 'y', 'size']
        new_list = []
        with open(table, 'r') as f:
            i = f.readlines()
            n = 0
            for line in i:
                n += 1
                if n == 1:
                    pass
                else:
                    list1 = line.strip().split('\t')
                    new_list.append(list1[0])
        new_df.insert(0, "name", new_list)
        new_df.to_csv(self.work_dir + "/new_table.xls", header=1, index=None, sep='\t')

    def set_output(self):
        """
        设置输出文件路径
        """
        path1 = self.output_dir + "/bubble.xls"
        if os.path.exists(path1):
            os.remove(path1)
        os.link(self.work_dir + "/new_table.xls", self.output_dir + "/bubble.xls")

    def set_db(self):
        self.logger.info("开始导表")
        api_bubble = self.api.api("tool_lab.bubble")
        self.logger.info('------------------')
        self.logger.info(self.output_dir + '/bubble/bubble.xls')
        api_bubble.add_bubble_detail(self.option('main_id'), self.output_dir + '/bubble/bubble.xls')
        self.logger.info("导表结束")
        self.end()
