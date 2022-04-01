# -*- coding: utf-8 -*-
# __author__ = 'linna'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError


class Regression1Agent(Agent):
    def __init__(self, parent):
        super(Regression1Agent, self).__init__(parent)
        options = [
            {"name": "tooltable", "type": "infile", "format": "toolapps.table"},
            {"name": "group_table", "type": "infile", "format": "toolapps.group_table"},
            {"name": "scale", "type": "string", "default": "F"},
            {"name": "themename", "type": "string", "default": "线性回归分析"},
            {"name": "main_id", "type": "string"},
            {"name": "specimen_name", "type": "string", "default": "列标签"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('tooltable').is_set:
            raise OptionError('必须提供数据表')
        self.option('tooltable').get_info()
        self.logger.info(self.option('specimen_name'))
        if self.option('group_table').is_set:
            if self.option('specimen_name') == '列标签':
                for i in self.option('group_table').prop['sample_name']:
                    if i not in self.option('tooltable').prop['col_sample']:
                        raise OptionError('分组文件中的样本不存在于表格中，查看是否数据选择错误', code="32702909")
            else:
                for i in self.option('group_table').prop['sample_name']:
                    if i not in self.option('tooltable').prop['row_sample']:
                        raise OptionError('分组文件中的样本不存在于表格中，查看是否数据选择错误', code="32702910")
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        if self.option('group_table'):
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "线性回归分析结果输出目录"],
                ["./regression_data.xls", "xls", "数据表"],
                ["./regression_message.xls", "xls", "信息表"],
                ["./group.xls", "xls", "分组表"],
            ])
        else:
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "线性回归分析结果输出目录"],
                ["./regression_data.xls", "xls", "数据表"],
                ["./regression_message.xls", "xls", "信息表"],
            ])
        super(Regression1Agent, self).end()


class Regression1Tool(Tool):
    def __init__(self, config):
        super(Regression1Tool, self).__init__(config)
        self.perl_path = "program/perl-5.24.0/bin/perl"
        self.script_path = self.config.PACKAGE_DIR + "/statistical/regression_analysis1.pl"
        self.R_path = "program/R-3.3.1/bin/Rscript"

    def run_regression(self):
        cmd = self.perl_path +' %s -i %s -scale %s -o %s ' % (self.script_path, self.formattable(self.option('tooltable').path), self.option('scale'), self.work_dir + '/regression')
        self.logger.info('运行regression_analysis1程序')
        self.logger.info(cmd)
        command = self.add_command('cmd', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("Regression的r文件生成成功")
        else:
            self.set_error("Regression的r文件生成失败")
            raise Exception("Regression的r文件生成失败")
        cmd_1 = '%s %s' % (self.R_path, self.work_dir + "/cmd.r")
        self.logger.info(cmd_1)
        command1 = self.add_command('cmd_1', cmd_1)
        command1.run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("R程序计算Regression成功")
        else:
            self.set_error("R程序计算Regression失败")
            raise Exception("R程序计算Regression失败")

    def run(self):
        super(Regression1Tool, self).run()
        self.run_regression()
        self.linkfile()
        self.set_db()
        self.end()

    def linkfile(self):
        if self.option('group_table').is_set:
            group_file = self.option('group_table').path
            group_link = self.output_dir + '/group.xls'
            if os.path.exists(group_link):
                os.remove(group_link)
            os.link(group_file, group_link)
            linksites_data = os.path.join(self.output_dir, 'regression_data.xls')
            linksites_r = os.path.join(self.output_dir, 'regression_message.xls')
            for i in linksites_data, linksites_r:
                if os.path.exists(i):
                    os.remove(i)
            os.link(self.work_dir + '/regression_data.xls', linksites_data)
            os.link(self.work_dir + '/regression_message.xls', linksites_r)
        else:
            linksites_data = os.path.join(self.output_dir, 'regression_data.xls')
            linksites_r = os.path.join(self.output_dir, 'regression_message.xls')
            for i in linksites_data, linksites_r:
                if os.path.exists(i):
                    os.remove(i)
            os.link(self.work_dir + '/regression_data.xls', linksites_data)
            os.link(self.work_dir + '/regression_message.xls', linksites_r)

    def formattable(self, tablepath):
        """
        转置表格
        """
        this_table = tablepath
        if self.option('specimen_name') != '列标签':
            return this_table
        else:
            newtable = this_table + '.T'
            self.t_table(this_table, newtable)
            return newtable

    def t_table(self, table_file, new_table):
        """
        转置表格实现
        """
        with open(table_file) as f, open(new_table, 'w') as w:
            table_list = [i.rstrip().split('\t') for i in f.readlines()]
            table_list = map(lambda *a: '\t'.join(a) + '\n', *table_list)
            w.writelines(table_list)

    def set_db(self):
        api_regression = self.api.api('tool_lab.regression')
        api_regression.add_regression(self.option('main_id'), self.output_dir, self.formattable(self.option('tooltable').path))
        self.logger.info("导表结束")


