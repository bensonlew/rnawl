# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
import os
import re


class Fasta2sqnAgent(Agent):
    def __init__(self, parent):
        super(Fasta2sqnAgent, self).__init__(parent)
        options = [
            {'name': 'temp_file', 'type': 'string'},
            {'name': 'seq_file', 'type': 'infile', 'format': 'sequence.fasta'},
            {'name': 'structure_type', 'type': 'string'},
            {'name': 'structure_file', 'type': 'string'},
            {'name': 'annot_file', 'type': 'infile', 'format':'sequence.profile_table'},
            {'name': 'species', 'type': 'string'},
            {'name': 'lab_name', 'type': 'string'},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option('temp_file'):
            self.set_error('缺少template文件！')
        if not self.option('seq_file').is_set:
            self.set_error('缺少输入fasta序列文件！')
        if not (self.option('structure_type') and self.option('structure_file')):
            self.set_error('必须同时设置结构文件和文件类型')

    def set_resource(self):
        self._cpu = 1
        self._memory = '30G'

    def end(self):
        super(Fasta2sqnAgent, self).end()


class Fasta2sqnTool(Tool):
    def __init__(self, config):
        super(Fasta2sqnTool, self).__init__(config)
        self.tbl2asn = self.config.SOFTWARE_DIR + '/bioinfo/tool_lab/tbl2asn/tbl2asn'
        #self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/library/glibc-2.14/lib/')
        self.env = 'export LD_LIBRARY_PATH={}:$LD_LIBRARY_PATH'.format(self.config.SOFTWARE_DIR + '/library/glibc-2.14/lib/')
        self.fasta2sqn = self.config.PACKAGE_DIR + '/tool_lab/fasta2sqn.py'
        self.python = '/miniconda2/bin/python'

    def run(self):
        super(Fasta2sqnTool, self).run()
        self.f2s_run()
        self.set_output()
        self.end()

    def f2s_run(self):
        cmd = '{} {} -s {} -f {}'.format(self.python, self.fasta2sqn,
                                         self.option('temp_file'), self.option('seq_file').prop['path'])
        if self.option('structure_type') == 'tbl':
            cmd += ' -t {} --tbl2asn "{}"'
        elif self.option('structure_type') == 'gtf':
            cmd += ' -g {} --tbl2asn "{}"'
        elif self.option('structure_type') == 'info':
            cmd += ' -i {} --tbl2asn "{}"'
        cmd = cmd.format(self.option('structure_file'), self.env + ' && ' + self.tbl2asn)
        cmd += ' -o ' + self.output_dir
        if self.option('annot_file').is_set:
            cmd += ' -a ' + self.option('annot_file').path
        if self.option('species'):
            cmd += ' -c ' + self.option('species')
        if self.option('lab_name'):
            cmd += ' -l ' + self.option('lab_name')
        command = self.add_command('run_fasta2sqn', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('fasta to sqn file done')
        else:
            self.logger.info('wrong')

    def set_output(self):
        self.logger.info('设置结果文件')
