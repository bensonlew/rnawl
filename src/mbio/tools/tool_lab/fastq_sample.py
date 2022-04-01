# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
import os
import re


class FastqSampleAgent(Agent):
    def __init__(self, parent):
        super(FastqSampleAgent, self).__init__(parent)
        options = [
            {'name': 'read1', 'type': 'infile', 'format': 'sequence.fastq'},
            {'name': 'read2', 'type': 'infile', 'format': 'sequence.fastq'},
            {'name': 'read', 'type': 'infile', 'format': 'sequence.fastq'},
            {'name': 'extract_type', 'type': 'string', },
            {'name': 'threshold', 'type': 'float', },
            {'name': 'nthread', 'type': 'int', 'default': 2},
        ]
        self.add_option(options)

    def check_options(self):
        if not (self.option('extract_type') or self.option('threshold')):
            raise OptionError('请设置收取条件')
        if self.option('extract_type') not in ['size', 'base', 'seq_num']:
            raise OptionError('目前可选的抽取方式只包括 size, base, seq')
        return True

    def set_resource(self):
        self._cpu = self.option('nthread')
        self._memory = '2G'


class FastqSampleTool(Tool):
    def __init__(self, config):
        super(FastqSampleTool, self).__init__(config)
        self.seqkit = '/bioinfo/seq/seqkit'

    def run(self):
        super(FastqSampleTool, self).run()
        self.run_sample()
        self.set_output()
        self.end()

    def run_sample(self):
        proportion = self.get_threshold()
        cmd = self.seqkit + ' sample -p {} -j {} '.format(proportion, self.option('nthread')) +\
            '-o {} {}'
        commands = []
        if self.option('read').is_set:
            prefix = os.path.basename(self.option('read').prop['path'])
            prefix = re.sub(r'\.(\w+)$', r'_ext.\1', prefix)
            cmd = cmd.format(self.output_dir + '/' + prefix, self.option('read').prop['path'])
            command = self.add_command('single_ext', cmd1).run()
            commands.append(command)
        else:
            prefix1 = os.path.basename(self.option('read1').prop['path'])
            prefix1 = re.sub(r'\.(\w+)$', r'_ext.\1', prefix1)
            prefix2 = os.path.basename(self.option('read2').prop['path'])
            prefix2 = re.sub(r'\.(\w+)$', r'_ext.\1', prefix2)
            cmd1 = cmd.format(self.output_dir + '/' + prefix1, self.option('read1').prop['path'])
            cmd2 = cmd.format(self.output_dir + '/' + prefix2, self.option('read2').prop['path'])
            command1 = self.add_command('read1_ext', cmd1).run()
            command2 = self.add_command('read2_ext', cmd2).run()
            commands.extend([command1, command2])
        self.wait()
        for c in commands:
            if c.return_code == 0:
                self.logger.info('{} done'.format(c.name))
            else:
                self.set_error('wrong in {}'.format(c.name))

    def get_threshold(self):
        if self.option('read').is_set:
            seq = self.option('read').prop['path']
        else:
            seq = self.option('read1').prop['path']
        if self.option('extract_type') == 'size':
            file_size = float(os.path.getsize(seq)) / (1024 * 1024 * 1024)
            prop = self.option('threshold') / file_size
        else:
            cmd = self.seqkit + ' stat {} -o tmp.stats'.format(seq)
            command = self.add_command('fq_stat', cmd).run()
            self.wait(command)
            if command.return_code == 0:
                with open('tmp.stats', 'r') as stat:
                    line = stat.readlines()[1].split()
                if self.option('extract_type') == 'base':
                    prop = self.option('threshold') / float(re.sub(r'\,', '', line[4]))
                else:
                    prop = self.option('threshold') / float(re.sub(r'\,', '', line[3]))
        if prop > 1:
            prop = 1
        return prop

    def set_output(self):
        pass
