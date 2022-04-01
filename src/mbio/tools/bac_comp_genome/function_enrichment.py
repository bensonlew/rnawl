# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import os
import pandas as pd

from biocluster.agent import Agent
from biocluster.config import Config
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool


class FunctionEnrichmentAgent(Agent):
    def __init__(self, parent):
        super(FunctionEnrichmentAgent, self).__init__(parent)
        options = [
            {'name': 'annotable', 'type': 'infile', 'format': 'sequence.profile_table'},
            {'name': 'levelenriched', 'type': 'string'},
            {'name': 'samples', 'type': 'string'},
            {'name': 'testset', 'type': 'string', 'default': ''},
            {'name': 'enrichout', 'type': 'string', 'default': ''},
            {'name': 'corrected', 'type': 'string'},
            {'name': 'go', 'type': 'string'},
            {'name': 'functype', 'type': 'string', 'default': 'KEGG'},
            {'name': 'formated', 'type': 'string', 'default': 'N'}
        ]
        self.add_option(options)
        self.step.add_steps('func_enrich')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.func_enrich.start()
        self.step.update()

    def stepfinish(self):
        self.step.func_enrich.finish()
        self.step.update()

    def check_options(self):
        pass

    def set_resource(self):
        self._memory = '4G'
        self._cpu = 2

    def end(self):

        super(FunctionEnrichmentAgent, self).end()


class FunctionEnrichmentTool(Tool):
    def __init__(self, config):
        super(FunctionEnrichmentTool, self).__init__(config)
        self.python = '/program/Python/bin/python'
        self.package_path = self.config.PACKAGE_DIR + '/bac_comp_genome/enrichment.py'

    def run(self):
        super(FunctionEnrichmentTool, self).run()
        self.run_enrich()
        self.set_output()
        self.end()

    def run_enrich(self):
        cmd = self.python + ' ' + self.package_path +\
            ' -a {} -s {} -o {} -c {} -l {} -g {} -t {} -f {}'.format(
                self.option('annotable').prop['path'], self.option('testset'),
                self.option('enrichout'), self.option('corrected'),
                self.option('levelenriched'), self.option('go'),
                self.option('functype'), self.option('formated')
            )
        if self.option('samples'):
            cmd += ' -sp ' + self.option('samples')
        command = self.add_command('function_enrich', cmd).run()
        self.wait(command)

        if command.return_code == 0:
            pass
        else:
            self.set_error('富集分析失败')

    def set_output(self):
        for root, dirs, files in os.walk(self.output_dir):
            for f in files:
                os.remove(os.path.join(root, f))
        self.logger.info('开始设置富集分析结果目录')
        output = '/' + self.option('enrichout')
        try:
            os.link(self.work_dir + output, self.output_dir + output)
            self.logger.info('成功设置富集分析结果文件')
        except:
            self.logger.info('设置富集分析结果文件失败')

