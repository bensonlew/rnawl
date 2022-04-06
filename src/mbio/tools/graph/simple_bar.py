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


class SimpleBarAgent(Agent):
    """
    version 1.0
    author: wangzhaoyue
    last_modify: 2017.04.25
    """
    def __init__(self, parent):
        super(SimpleBarAgent, self).__init__(parent)
        options = [
            {"name": "input_table", "type": "infile", "format": "toolapps.table"},  # 输入的表格，矩阵
            {"name": "group_table", "type": "infile", "format": "toolapps.group_table"},  # 输入的group表格
            {"name": "method", "type": "string", "default": "row"},  # 样本名的方向，默认样本在行row,column
            {"name": "combined_value", "type": "string", "default": "0.01"},  # 合并小于此值的属性
            {"name": "calculation", "type": "string", "default": "none"}  # 组内合并参数，none,sum,average,middle
        ]
        self.add_option(options)
        self.step.add_steps('simple_bar')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.simple_bar.start()
        self.step.update()

    def step_end(self):
        self.step.simple_bar.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        :return:
        """
        if not self.option("input_table"):
            raise OptionError("参数input_table不能为空")
        if self.option('group_table').is_set:
            if self.option('method') == 'column':
                for i in self.option('group_table').prop['sample_name']:
                    if i not in self.option('input_table').prop['row_sample']:
                        raise Exception('分组文件中的样本{}不存在于表格第一列中，查看是否是数据取值选择错误'.format(i))
            else:
                for i in self.option('group_table').prop['sample_name']:
                    if i not in self.option('input_table').prop['col_sample']:
                        raise Exception('分组文件中的样本{}不存在于表格第一行中，查看是否是数据取值选择错误'.format(i))

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 15
        self._memory = '1G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "柱形图结果目录"],
            # ["./.*final_value.xls", "xls", "结果表"],
        ])
        super(SimpleBarAgent, self).end()


class SimpleBarTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(SimpleBarTool, self).__init__(config)
        self._version = 1.0
        self.Python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/python '
        self.Python_path2 = '/miniconda2/bin/python '
        self.path = self.config.PACKAGE_DIR + '/statistical/simple_operation.py'
        self.software = 'program/parafly-r2013-01-21/bin/bin/ParaFly'

    def create_common_table(self):
        """
        输入的文件统一处理成标准格式的文件,第一列为样本名
        """
        self.logger.info("开始分析")
        cmd_list = []
        input_table = self.option("input_table").prop['new_table']
        # 有分组方案时，判断组内合并方式，建立每个分组方案的命令
        if self.option("group_table").is_set and self.option("calculation") != "none":
            for i in range(len(self.option("group_table").prop['group_scheme'])):
                select_group = []
                sample_dir = self.work_dir + '/' + self.option("group_table").prop['group_scheme'][i]
                os.mkdir(sample_dir)
                select_group.append(self.option("group_table").prop['group_scheme'][i])
                self.option('group_table').sub_group(sample_dir + '/group_' + str(i + 1), select_group)  # 获得用于计算的一个分组方案
                # group_table = sample_dir + '/group_' + str(i + 1)
                os.link(input_table, sample_dir + '/input_' + str(i + 1))
                # self.option("input_table").get_table_of_main_table(input_table, sample_dir + '/input_' + str(i + 1),
                #                                                    group_table)
                new_input_table = sample_dir + '/input_' + str(i + 1)
                middle_input = sample_dir + "/" + self.option("group_table").prop['group_scheme'][
                    i] + '_middle_input.xls'
                final_input = sample_dir + "/" + self.option("group_table").prop['group_scheme'][
                    i] + "_final_input.xls"  # 样本在列，方便计算
                combined_txt = sample_dir + "/" + self.option("group_table").prop['group_scheme'][
                    i] + "_final_table.xls"
                value_table = sample_dir + "/" + self.option("group_table").prop['group_scheme'][
                    i] + "_final_value.xls"
                cmd = self.Python_path + self.path + " -i %s -method %s -i1 %s -i2 %s -o1 %s -o2 %s -combined %s -group %s -calculation %s" % (
                    new_input_table, self.option('method'), middle_input, final_input, combined_txt, value_table,
                    self.option("combined_value"),
                    sample_dir + '/group_' + str(i + 1), self.option("calculation"))
                self.logger.info(cmd)
                cmd_list.append(cmd)
            # 循环投递
            self.logger.info('运行python脚本，进行计算')
            n = len(cmd_list) / 15
            if len(cmd_list) % 15 != 0:
                n += 1
            for i in range(0, n):
                cmd_file = os.path.join(self.work_dir + '/' + self.option("group_table").prop['group_scheme'][i],
                                        'cmd_list_{}.txt'.format(i + 1))
                wrong_cmd = os.path.join(self.work_dir + '/' + self.option("group_table").prop['group_scheme'][i],
                                         'failed_cmd_{}.txt'.format(i + 1))
                with open(cmd_file, "w")as c:
                    for n in range(0, 15):
                        if len(cmd_list) == 0:
                            break
                        cmd = cmd_list.pop()
                        c.write(cmd + "\n")
                final_cmd = '{} -c {} -CPU 10 -failed_cmds {}'.format(self.software, cmd_file, wrong_cmd)
                command = self.add_command("cmd_{}".format(i + 1), final_cmd).run()
                self.wait(command)
                if command.return_code == 0:
                    self.logger.info("运行{}完成".format(command.name))
                else:
                    self.set_error("运行{}运行出错!".format(command.name))
                    raise Exception("运行柱形图运行出错，请检查输入的表格是否正确")
        else:
            middle_input = self.work_dir + "/middle_input.xls"  # 样本在行
            final_input = self.work_dir + "/final_input.xls"  # 样本在列，方便计算
            combined_txt = self.work_dir + "/final_table.xls"
            value_table = self.work_dir + "/final_value.xls"
            cmd = self.Python_path2 + self.path + " -i %s -method %s -i1 %s -i2 %s -o1 %s -o2 %s -combined %s" % (
                input_table, self.option('method'), middle_input, final_input, combined_txt, value_table, self.option("combined_value"))
            command = self.add_command("cmd_bar", cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("运行{}完成".format(command.name))
            else:
                self.set_error("运行{}运行出错!".format(command.name))
                raise Exception("运行柱形图运行出错，请检查输入的表格是否正确")

    def set_output(self):
        """
        将结果文件链接至output
        """
        self.logger.info("set output")
        if self.option("group_table").is_set and self.option("calculation") != "none":
            for i in range(len(self.option("group_table").prop['group_scheme'])):
                input_path = self.work_dir + '/' + self.option("group_table").prop['group_scheme'][i] + "/" + \
                             self.option("group_table").prop['group_scheme'][i]
                output_path = self.output_dir + "/" + self.option("group_table").prop['group_scheme'][i]
                os.link(input_path + '_final_value.xls', output_path + '_final_value.xls')  # 丰度表格
                os.link(input_path + '_final_table.xls', output_path + '_matrix_bar.xls')  # 百分比表格
        else:
            os.link(self.work_dir + '/final_value.xls', self.output_dir + '/final_value.xls')  # 丰度表格
            os.link(self.work_dir + '/final_table.xls', self.output_dir + '/matrix_bar.xls')  # 百分比表格
        self.logger.info("done")

    def run(self):
        """
        运行
        """
        super(SimpleBarTool, self).run()
        self.create_common_table()
        self.set_output()
        self.end()
