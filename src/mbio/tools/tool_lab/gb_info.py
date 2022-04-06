# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
import os


class GbInfoAgent(Agent):
    def __init__(self, parent):
        super(GbInfoAgent, self).__init__(parent)
        options = [
            {'name': 'gb_file', 'type': 'infile', 'format': 'gene_structure.gbk'},
            {'name': 'cds', 'type': 'int', 'default': 0}, 
            {'name': 'genome', 'type': 'int', 'default': 0},
            {'name': 'prot', 'type': 'int', 'default': 0},
            {'name': 'gff', 'type': 'int', 'default': 0},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option('gb_file').is_set:
            raise OptionError('请正确输入文本文件')
        if sum([self.option('gff'), self.option('prot'), self.option('cds'), self.option('genome')]) == 0:
            raise OptionError('请至少选择一种需要输出的格式: cds, prot, gff, genome')
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '2G'

    def end(self):
        super(GbInfoAgent, self).end()



class GbInfoTool(Tool):
    def __init__(self, config):
        super(GbInfoTool, self).__init__(config)
        self.gb_info = self.config.PACKAGE_DIR + '/tool_lab/gb_info.py'
        self.python = '/miniconda2/bin/python'

    def run(self):
        super(GbInfoTool, self).run()
        self.run_gb_info()
        self.set_output()
        self.end()

    def run_gb_info(self):
        cmd = self.python + ' {0} {1}'
        cmd = cmd.format(self.gb_info, self.option('gb_file').prop['path'])
        if self.option('gff'):
            cmd += ' -gff'
        if self.option('cds'):
            cmd += ' -cds'
        if self.option('prot'):
            cmd += ' -prot'
        if self.option('genome'):
            cmd += ' -genome'
        command = self.add_command('gb_info', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info('gb_info done')
        else:
            self.set_error('wrong in gb_info')

    def set_output(self):
        prefix = os.path.splitext(os.path.basename(self.option('gb_file').path))[0]
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        for f in os.listdir(self.work_dir):
            if f.startswith(prefix):
                os.link(f, os.path.join(self.output_dir, f))

