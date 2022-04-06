# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'
import os
import shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.core.exceptions import FileError
import re
from collections import defaultdict


class BoxPlotAgent(Agent):
    """
    version 1.0
    author: wangzhaoyue
    last_modify: 2017.04.25
    """
    def __init__(self, parent):
        super(BoxPlotAgent, self).__init__(parent)
        options = [
            {"name": "input_table", "type": "infile", "format": "toolapps.table"},  # 输入的表格，矩阵
            {"name": "method", "type": "string", "default": "row"},  # 统计数据的方向，row,column
            {"name": "group_table", "type": "infile", "format": "toolapps.group_table"},  # 输入的分组表格，矩阵
            {"name": "first_group", "type": "string", "default": ""},  # 一级分组名称
            {"name": "sed_group", "type": "string", "default": ""},  # 二级分组名称
        ]
        self.add_option(options)
        self.step.add_steps('box_plot')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.box_plot.start()
        self.step.update()

    def step_end(self):
        self.step.box_plot.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        :return:
        """
        if not self.option("input_table"):
            raise OptionError("参数input_table不能为空")
        if self.option('group_table').is_set and self.option("sed_group") != '':
            if self.option("first_group") == '':
                raise OptionError("必须填写一级分组方案名称，才能填写二级分组方案名称")
            else:
                pass
            if self.option('first_group') not in self.option("group_table").prop['group_scheme']:
                raise OptionError("填写的分组方案名称%s在分组文件中不存在，请核实！"%(self.option('first_group')))
            if self.option('sed_group') not in self.option("group_table").prop['group_scheme']:
                raise OptionError("填写的分组方案名称%s在分组文件中不存在，请核实！" % (self.option('sed_group')))
        if self.option('group_table').is_set:
            if self.option('method') == 'column':
                for i in self.option('group_table').prop['sample_name']:
                    if i not in self.option('input_table').prop['row_sample']:
                        raise OptionError('分组文件中的样本{}不存在于表格第一列中，查看是否是数据取值选择错误'.format(i))
            else:
                for i in self.option('group_table').prop['sample_name']:
                    if i not in self.option('input_table').prop['col_sample']:
                        raise OptionError('分组文件中的样本{}不存在于表格第一行中，查看是否是数据取值选择错误'.format(i))

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "箱线图结果目录"],
        ])
        super(BoxPlotAgent, self).end()


class BoxPlotTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(BoxPlotTool, self).__init__(config)
        self._version = 1.0
        self.Python_path = '/miniconda2/bin/python '
        self.box_path = self.config.SOFTWARE_DIR + '/bioinfo/statistical/scripts/boxplot.py '

    def create_common_table(self):
        """
        输入的文件统一处理成标准格式的文件,第一列为样本名
        """
        final_input = self.work_dir + '/input_table.xls'
        with open(self.option("input_table").prop["new_table"]) as f, open(final_input, "w+")as fw:
            lines = f.readlines()
            if self.option("method") == "column":
                table_list = [i.rstrip().split('\t') for i in lines]
                lines = map(lambda *a: '\t'.join(a) + '\n', *table_list)
            for i in lines:
                self.logger.info(i)
                fw.write(i)

    def box_plot_run(self):
        cmd = self.Python_path + self.box_path + " %s" % (self.work_dir + '/input_table.xls')
        command = self.add_command("cmd_box", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行{}完成".format(command.name))
        else:
            self.set_error("运行{}运行出错!".format(command.name))
            raise Exception("运行箱线图运行出错，请检查输入的表格是否正确")
        self.set_output()

    def set_output(self):
        """
        将结果文件链接至output
        """
        self.logger.info("set output")
        shutil.copy2(self.work_dir + "/test.boxplot.xls", self.output_dir + "/boxplot.xls")
        self.logger.info("done")

    def run(self):
        """
        运行
        """
        super(BoxPlotTool, self).run()
        if self.option('group_table').is_set:
            if self.option("sed_group") == '':
                group_list = self.option("group_table").prop['group_scheme']
                for i in range(len(group_list)):
                    select_list = []
                    select_list.append(group_list[i])
                    self.option('group_table').sub_group(self.output_dir + '/group_' + str(i + 1), select_list)
            else:
                first_list = []
                second_list = []
                first_list.append(self.option("first_group"))
                second_list.append(self.option("sed_group"))
                self.option('group_table').sub_group(self.output_dir + '/first_group.xls', first_list)
                self.option('group_table').sub_group(self.output_dir + '/sed_group.xls', second_list)
        self.create_common_table()
        self.box_plot_run()
        self.end()
