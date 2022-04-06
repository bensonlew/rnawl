# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
import os
import pandas as pd

from biocluster.agent import Agent
from biocluster.config import Config
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool


class KeggGraphInfoAgent(Agent):
    def __init__(self, parent):
        super(KeggGraphInfoAgent, self).__init__(parent)
        options = [
            {'name': 'annotable', 'type': 'infile', 'format': 'sequence.profile_table'},
            # {'name': 'graph_info', 'type': 'string'},
            # {'name': 'out', 'type': 'string'},
        ]
        self.add_option(options)
        self.step.add_steps('kegg_graph')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.kegg_graph.start()
        self.step.update()

    def stepfinish(self):
        self.step.kegg_graph.finish()
        self.step.update()

    def check_options(self):
        if not self.option('annotable').is_set:
            raise OptionError('必须设置输入文件')
        return True

    def set_resource(self):
        self._memory = '4G'
        self._cpu = 2

    def end(self):

        super(KeggGraphInfoAgent, self).end()


class KeggGraphInfoTool(Tool):
    def __init__(self, config):
        super(KeggGraphInfoTool, self).__init__(config)
        self.python = '/miniconda2/bin/python'
        self.package_path = self.config.PACKAGE_DIR +\
            '/bac_comp_genome/kegg_graph_info.py'
        self._graph_info = self.config.SOFTWARE_DIR +\
            '/database/KEGG/bac_comp/kegg_KOinMap_info.xls'

    def run(self):
        super(KeggGraphInfoTool, self).run()
        self.run_graph()
        self.set_output()
        self.end()

    def run_graph(self):
        cmd = self.python + ' ' + self.package_path +\
            ' -a {} -g {}'.format(self.option('annotable').prop['path'], self._graph_info)
        command = self.add_command('kegg_graph', cmd).run()
        self.wait(command)

        if command.return_code == 0:
            pass
        else:
            self.set_error('kegg graph运行失败')

    def set_output(self):
        for root, dirs, files in os.walk(self.output_dir):
            for f in files:
                os.remove(os.path.join(root, f))
        self.logger.info('开始设置kegg graph结果目录')
        output = '/kegg_graph_info.xls'
        try:
            os.link(self.work_dir + output, self.output_dir + output)
            self.logger.info('成功设置kegg graph结果文件')
        except:
            self.logger.info('设置kegg graph结果文件失败')

