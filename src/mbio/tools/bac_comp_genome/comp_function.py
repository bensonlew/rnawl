# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import os
import pandas as pd

from biocluster.agent import Agent
from biocluster.config import Config
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool


class CompFunctionAgent(Agent):
    def __init__(self, parent):
        super(CompFunctionAgent, self).__init__(parent)
        options = [
            {'name': 'annotable', 'type': 'string'},
            {'name': 'functiontype', 'type': 'string'},
            {'name': 'level', 'type': 'string'},
            {'name': 'splist', 'type': 'string'},
            {'name': 'groups', 'type': 'string', 'default': 'ALL'},
            {'name': 'corepan', 'type': 'string', 'default': ''},
            {'name': 'pancat', 'type': 'string', 'default': ''},
            {'name': 'output', 'type': 'string', 'default': 'out'},
            {'name': 'selectedcat', 'type': 'string', 'default': 'ALL'},
            {'name': 'abundance', 'type': 'string', 'default': 'F'},
            {'name': 'corepansets', 'type': 'string', 'default': ''},
            {'name': 'result_type', 'type': 'string', 'default': ''}
        ]
        self.add_option(options)
        self.step.add_steps('comp_function')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.comp_function.start()
        self.step.update()

    def stepfinish(self):
        self.step.comp_function.finish()
        self.step.update()

    def check_options(self):
        if not self.option('corepansets'):
            if not self.option('annotable'):
                raise OptionError('请输入注释表')
            if not self.option('functiontype'):
                raise OptionError('必须指定比较的功能类别')
            if not self.option('splist'):
                raise OptionError('请给出样本分组文件')

    def set_resource(self):
        self._cpu = 2
        self._memory = '4G'

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules(['.', "", "结果输出目录"])
        super(CompFunctionAgent, self).end()


class CompFunctionTool(Tool):
    def __init__(self, config):
        super(CompFunctionTool, self).__init__(config)
        #self.python = self.config.SOFTWARE_DIR + '/program/Python/bin/python'
        self.python = '/program/Python/bin/python'
        self.package_path = self.config.PACKAGE_DIR + '/bac_comp_genome/comp_function.py'

    def run(self):
        super(CompFunctionTool, self).run()
        self.comp_run()
        self.set_output()
        self.end()

    def comp_run(self):
        cmd = self.python + ' ' + self.package_path +\
            ' -a {} -f {} -l {} -s {} -o {} -b {} -g {}'.format(
                self.option('annotable'), self.option('functiontype'),
                self.option('level'), self.option('splist'),
                self.option('output'), self.option('abundance'),
                self.option('groups')
            )
        if self.option('corepan'):
            if self.option('corepansets'):
                cmd += " -c {} -p '{}' --corepansets {}".format(
                        self.option('corepan'), self.option('pancat'),
                        self.option('corepansets')
                        )
            else:
                cmd += " -c {} -p '{}' -d {}".format(
                        self.option('corepan'), self.option('pancat'),
                        self.option('selectedcat')
                    )
            if self.option('result_type'):
                cmd += " -r " + self.option('result_type')
        command = self.add_command('func_comp', cmd, ignore_error=True).run()
        self.wait(command)

        if command.return_code == 0:
            self.logger.info('package comp_funcion.py运行成功')
        else:
            with open('func_comp.o', 'r') as o:
                for l in o:
                    if l.startswith('Exception'):
                        self.set_error(l)
            self.set_error('package comp_funcion.py运行出错')

    def set_output(self):
        for root, dirs, files in os.walk(self.output_dir):
            for f in files:
                os.remove(os.path.join(root, f))
        self.logger.info('开始设置结果输出目录')
        output = '/' + self.option('output')
        cp = '/corepan_sets.xls'
        try:
            if self.option('corepansets'):
                os.link(self.work_dir + cp, self.output_dir + cp)
            else:
                os.link(self.work_dir + output, self.output_dir + output)
            self.logger.info('成功设置结果文件到结果目录')
        except:
            self.set_error('设置结果文件到目录失败')

