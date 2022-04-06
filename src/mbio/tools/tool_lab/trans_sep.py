# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
import os


class TransSepAgent(Agent):
    def __init__(self, parent):
        super(TransSepAgent, self).__init__(parent)
        options = [
            {'name': 'input', 'type': 'string',},
            {'name': 'output', 'type': 'string'},
            {'name': 'sep_from', 'type': 'string'},
            {'name': 'sep_to', 'type': 'string'},
        ]
        self.add_option(options)

    def check_options(self):
        if not os.path.isfile(self.option('input')):
            raise OptionError('请正确输入文本文件')
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '2G'


class TransSepTool(Tool):
    def __init__(self, config):
        super(TransSepTool, self).__init__(config)
        self.trans_sep = self.config.PACKAGE_DIR + '/tool_lab/change_sign.py'
        self.python = '/miniconda2/bin/python'

    def run(self):
        super(TransSepTool, self).run()
        self.run_trans()
        self.set_output()
        self.end()

    def run_trans(self):
        sper = {
            'space': ' ',
            'comma': ',',
            'semicolon': ';',
            'underline': '_',
            'tab': '\t',
        }
        if self.option('sep_from') not in sper or self.option('sep_to') not in sper:
            self.set_error('可选的分隔符标识符{}'.format(sper.keys()))
        cmd = self.python + ' {} -i {} -o {} --insig "{}" --outsig "{}"'
        cmd = cmd.format(self.trans_sep, self.option('input'),
                         os.path.join(self.output_dir, self.option('output')),
                         sper[self.option('sep_from')], sper[self.option('sep_to')])
        command = self.add_command('trans_sep', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('trans done')
        else:
            self.set_error('wrong in trans sep')

    def set_output(self):
        pass
