# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
import subprocess
from biocluster.core.exceptions import OptionError



class RegressionAgent(Agent):
    """
    version v1.0
    author: gaohao
    last_modified:2017.10.30
    """

    def __init__(self, parent):
        super(RegressionAgent, self).__init__(parent)
        options = [
            {"name": "abu_table", "type": "infile","format": "meta.otu.otu_table,toolapps.table_regression"},
            {"name":"group","type":"infile","format": "toolapps.group_table"},
            {"name":"num_method","type":"string","default":"none"},
            {"name": "abu_trans", "type": "string", "default": "column"}
        ]
        self.add_option(options)
        self.step.add_steps('regression')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.regression.start()
        self.step.update()

    def step_end(self):
        self.step.regression.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('abu_table').is_set:
            raise OptionError('必须提供丰度表')
        if not self.option('group').is_set:
            raise OptionError('必须提供分组表')
        else:
            sample_e =[]
            sample_o =[]
            with open(self.option('group').prop['path'], 'r') as e:
                for line in e:
                    line = line.strip('\n\r').split('\t')
                    sample_e.append(line[0])
            sample_e.pop(0)
            with open(self.option('abu_table').prop['path'], 'r') as o:
                for line in o:
                    line = line.strip('\n\r').split('\t')
                    sample_o.append(line[0])
            sample_o.pop(0)
            sample_e.sort()
            sample_o.sort()
            for i in sample_o:
                if i in sample_e:
                    continue
                else:
                    raise OptionError('数据表中的样本和分组表中的样本不一致')

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
        ])
        super(RegressionAgent, self).end()


class RegressionTool(Tool):
    def __init__(self, config):
        super(RegressionTool, self).__init__(config)
        self._version = '1.0'
        self.perl_path = "program/perl-5.24.0/bin/perl"
        self.script1_path = self.config.PACKAGE_DIR + "/statistical/regression_analysis.pl"
        self.R_path = 'program/R-3.3.1/bin/Rscript'

    def run_regression(self):
        if self.option('abu_trans') in ['row']:
            abu_table = self.option('abu_table').prop['path']+ '.new'
            self.t_table(self.option('abu_table').prop['path'],abu_table)
        else:
            abu_table =self.option('abu_table').prop['path']
        if self.option('num_method') in ['center','scale','none']:
            cmd = self.perl_path +' %s -i %s -method %s -o %s ' % (self.script1_path, abu_table, self.option("num_method"),self.work_dir + '/Regression')
        else:
            raise OptionError('必须选择正确的数据处理方法！')
        self.logger.info(cmd)
        command = self.add_command('cmd', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("Regression的r文件生成成功")
        else:
            self.set_error("Regression的r文件生成失败")
            raise Exception("Regression的r文件生成失败")
        cmd_1 = '%s %s' % (self.R_path, self.work_dir + "/RegressionAnalysis.cmd.r")
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
        super(RegressionTool, self).run()
        self.run_regression()
        self.set_output()

    def set_output(self):
        linksites_data = os.path.join(self.output_dir, 'Regression.data.xls')
        linksites_r = os.path.join(self.output_dir, 'Regression.message.xls')
        for i in linksites_data,linksites_r:
            if os.path.exists(i):
               os.remove(i)
        os.link(self.work_dir + '/Regression.data.xls', linksites_data)
        os.link(self.work_dir + '/Regression.message.xls', linksites_r)
        self.end()

    def t_table(self, table_file, new_table):  # 表格转置
        """
    	转换颠倒表格内容
    	"""
        with open(table_file) as f, open(new_table, 'w') as w:
            table_list = [i.rstrip().split('\t') for i in f.readlines()]
            table_list = map(lambda *a: '\t'.join(a) + '\n', *table_list)
            w.writelines(table_list)



