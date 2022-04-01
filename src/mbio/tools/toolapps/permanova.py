# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
import subprocess
from biocluster.core.exceptions import OptionError


class PermanovaAgent(Agent):
    """
    脚本permanova.pl
    version v1.0
    author: gaohao
    last_modified:2017.10.30
    """

    def __init__(self, parent):
        super(PermanovaAgent, self).__init__(parent)
        options = [
            {"name": "abu_table", "type": "infile","format": "meta.otu.otu_table, meta.otu.tax_summary_dir, toolapps.table"},
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},
            {"name":"dis_method","type":"string","default":"bray"},
            {"name":"replacement_num","type":"int","default":"999"},
            {"name": "abu_trans", "type": "string", "default": "column"},
            {"name": "env_trans", "type": "string", "default": "column"}
        ]
        self.add_option(options)
        self.step.add_steps('permanova')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.permanova.start()
        self.step.update()

    def step_end(self):
        self.step.permanova.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('abu_table').is_set:
            raise OptionError('必须提供数据列表')
        if not self.option('envtable').is_set:
            raise OptionError('必须提供因子表')

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
        super(PermanovaAgent, self).end()


class PermanovaTool(Tool):
    def __init__(self, config):
        super(PermanovaTool, self).__init__(config)
        self._version = '1.0'
        self.perl_path = "program/perl-5.24.0/bin/perl"
        self.script1_path = self.config.PACKAGE_DIR + "/beta_diversity/permanova.pl"
        self.R_path = 'program/R-3.3.1/bin/Rscript'

    def run_permanova(self):
        if self.option('abu_trans')  in ['row']:
            abutable = self.option('abu_table').prop['path']+'.new'
            self.t_table(self.option('abu_table').prop['path'],abutable)            
        else:
            abutable = self.option('abu_table').prop['path']
        cmd = self.perl_path
        if self.option('env_trans') in ['row']:
            envtable = self.option('envtable').prop['path']+'.new'
            self.t_table(self.option('envtable').prop['path'],envtable)
        else:
            envtable = self.option('envtable').prop['path']
        if '_binary' in self.option('dis_method'):
            meth = re.search(r'(.*)_binary',self.option('dis_method'),re.M)
            cmd += ' %s -i %s -env %s -o %s -pe %s -binary T' % (
            self.script1_path, abutable,envtable,self.work_dir + '/permanova.xls',meth.group(1))
        else:
            cmd += ' %s -i %s -env %s -o %s -pe %s' % (
            self.script1_path,abutable,envtable,self.work_dir + '/permanova.xls',self.option('dis_method'))
        self.logger.info(cmd)
        command = self.add_command('cmd', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("permanova的r文件生成成功")
        else:
            self.set_error("permanova的r文件生成失败")
            raise Exception("permanova的r文件生成失败")
        cmd_1 = '%s %s' % (self.R_path, self.work_dir + "/cmd.r")
        self.logger.info(cmd_1)
        command1 = self.add_command('cmd_1', cmd_1)
        command1.run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("R程序计算permanova成功")
        else:
            self.set_error("R程序计算permanova失败")
            raise Exception("R程序计算permanova失败")

    def run(self):
        super(PermanovaTool, self).run()
        self.run_permanova()
        self.set_output()

    def set_output(self):
        linksites = os.path.join(self.output_dir, 'permanova.xls')
        if os.path.exists(linksites):
            os.remove(linksites)
        os.link(self.work_dir + '/permanova.xls', linksites)
        self.end()

    def t_table(self, table_file, new_table):  # 表格转置
        """
    	转换颠倒表格内容
    	"""
        with open(table_file) as f, open(new_table, 'w') as w:
            table_list = [i.rstrip().split('\t') for i in f.readlines()]
            table_list = map(lambda *a: '\t'.join(a) + '\n', *table_list)
            w.writelines(table_list)



