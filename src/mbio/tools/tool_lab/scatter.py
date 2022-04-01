# -*- coding: utf-8 -*-
# __author__ = 'linna'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import pandas as pd


class ScatterAgent(Agent):
    def __init__(self, parent):
        super(ScatterAgent, self).__init__(parent)
        options = [
            {"name": "tooltable", "type": "infile", "format": "tool_lab.scatter_table"},
            {"name": "x_axis", "type": "string"},
            {"name": "y_axis", "type": "string"},
            {"name": "z_axis", "type": "string"},
            {"name": "group", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数检查
        """
        self.logger.info(self.option('tooltable'))
        if not self.option('tooltable').is_set:
            raise OptionError('必须提供数据表', code="32702903")
        self.option('tooltable').get_info()
        if self.option('tooltable').prop['sample_num'] < 2:
            raise OptionError('列数少于2，不可进行分析', code="32702904")
        if self.option('x_axis') not in self.option('tooltable').prop['col_sample']:
            raise OptionError('x轴数据不在输入表格中，查看是否数据选择错误', code="32702909")
        if self.option('y_axis') not in self.option('tooltable').prop['col_sample']:
            raise OptionError('y轴数据不在输入表格中，查看是否数据选择错误', code="32702909")
        if self.option('z_axis'):
            if self.option('z_axis') not in self.option('tooltable').prop['col_sample']:
                raise OptionError('z轴数据不在输入表格中，查看是否数据选择错误', code="32702909")
        if self.option('group'):
            if self.option('group') not in self.option('tooltable').prop['col_sample']:
                raise OptionError('分组数据不在输入表格中，查看是否数据选择错误', code="32702909")
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(ScatterAgent, self).end()


class ScatterTool(Tool):
    def __init__(self, config):
        super(ScatterTool, self).__init__(config)
        self.table = self.get_table()

    def run(self):
        """
        运行
        """
        super(ScatterTool, self).run()
        self.get_table()
        self.linkfile()
        # self.set_db()
        self.end()

    def get_table(self):
        result = pd.read_table(self.option('tooltable').path)
        col_list = result.columns

        if self.option('z_axis'):
            new_result = result[[col_list[0], self.option('x_axis'), self.option('y_axis'), self.option('z_axis')]]
            new_result.rename(columns={col_list[0]: 'name', self.option('x_axis'): 'x', self.option('y_axis'): 'y',
                                       self.option('z_axis'): 'z'}, inplace=True)
            new_result.to_csv(os.path.join(self.work_dir, 'scatter.xls'), sep='\t', index=0)
            if self.option('group'):
                result = pd.read_table(self.option('tooltable').path)
                new_result = result[[col_list[0], self.option('x_axis'), self.option('y_axis'), self.option('z_axis'),
                                     self.option('group')]]
                new_result.rename(columns={col_list[0]: 'name', self.option('x_axis'): 'x', self.option('y_axis'): 'y',
                                           self.option('z_axis'): 'z', self.option('group'): 'group'}, inplace=True)
                new_result.to_csv(os.path.join(self.work_dir, 'scatter.xls'), sep='\t', index=0)
        else:
            new_result = result[[col_list[0], self.option('x_axis'), self.option('y_axis')]]
            new_result.rename(columns={col_list[0]: 'name', self.option('x_axis'): 'x',
                                       self.option('y_axis'): 'y'}, inplace=True)
            new_result.to_csv(os.path.join(self.work_dir, 'scatter.xls'), sep='\t', index=0)
            if self.option('group'):
                result = pd.read_table(self.option('tooltable').path)
                new_result = result[[col_list[0], self.option('x_axis'), self.option('y_axis'), self.option('group')]]
                new_result.rename(columns={col_list[0]: 'name', self.option('x_axis'): 'x', self.option('y_axis'): 'y',
                                           self.option('group'): 'group'}, inplace=True)
                new_result.to_csv(os.path.join(self.work_dir, 'scatter.xls'), sep='\t', index=0)
        return self.work_dir + '/scatter.xls'

    def linkfile(self):
        t_file = self.table
        link = self.output_dir + '/scatter.xls'
        if os.path.exists(link):
            os.remove(link)
        os.link(t_file, link)

    def set_db(self):
        self.logger.info("开始导表")
        api_scatter = self.api.api("tool_lab.scatter_api")
        api_scatter.add_scatter_detail(self.option('main_id'), self.output_dir)
        self.logger.info("导表结束")
