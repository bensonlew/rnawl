# -*- coding: utf-8 -*-
# __author__ = 'qindanhua,qinjincheng'

from __future__ import division

import os

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool

from mbio.packages.medical_transcriptome.functions import runcmd


class FastqDupAgent(Agent):
    def __init__(self, parent):
        super(FastqDupAgent, self).__init__(parent)
        options = [
            {'name': 'fastq_dir', 'type': 'infile', 'format': 'sequence.fastq_dir'},  # fastq文件夹
            {'name': 'fastq_s', 'type': 'infile', 'format': 'sequence.fastq'},  # 输入文件SE序列
            {'name': 'fastq_r', 'type': 'infile', 'format': 'sequence.fastq'},  # 输入文件PE的右端序列
            {'name': 'fastq_l', 'type': 'infile', 'format': 'sequence.fastq'},  # PE的左端序列
            {'name': 'fq_type', 'type': 'string'},
            {'name': 'sample_name', 'type': 'string', 'default': None},
            {'name': 'result', 'type': 'outfile', 'format': 'ref_rna_v2.common'}
        ]
        self.add_option(options)
        self._memory_increase_step = 150

    def check_options(self):
        if not self.option('fastq_dir').is_set and not self.option('fastq_r').is_set and not self.option(
                'fastq_s').is_set:
            raise OptionError('请传入fastq序列文件或者文件夹', code='34000901')
        if self.option('fq_type') not in ['PE', 'SE']:
            raise OptionError('请说明序列类型，PE or SE?', code='34000902')
        if not self.option('fastq_dir').is_set and self.option('fq_type') in ['PE']:
            if not self.option('fastq_r').is_set:
                raise OptionError('请传入PE右端序列文件', code='34000903')
            if not self.option('fastq_l').is_set:
                raise OptionError('请传入PE左端序列文件', code='34000904')
        if not self.option('fastq_dir').is_set:
            if self.option('fq_type') in ['SE'] and not self.option('fastq_s').is_set:
                raise OptionError('请传入SE序列文件', code='34000905')

    def set_resource(self):
        self._cpu = 1
        self._memory = '80G'

    def end(self):
        super(FastqDupAgent, self).end()


class FastqDupTool(Tool):
    def __init__(self, config):
        super(FastqDupTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python'
        }
        self.script = {
            'fastq_dup': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/seq/scripts/fastq_dup.py')
        }
        self.file = {
            'result': os.path.join(self.output_dir, '{}.fastq_dup.xls'.format(self.option('sample_name')))
        }

    def run(self):
        super(FastqDupTool, self).run()
        self.run_fastq_dup()
        self.set_output()
        self.end()

    def run_fastq_dup(self):
        cmd = '{} {}'.format(self.program['python'], self.script['fastq_dup'])
        if self.option('fq_type') == 'PE':
            cmd += ' -l {}'.format(self.option('fastq_l').prop['path'])
            cmd += ' -r {}'.format(self.option('fastq_r').prop['path'])
        elif self.option('fq_type') == 'SE':
            cmd += ' -s {}'.format(self.option('fastq_s').prop['path'])
        cmd += ' -o {}'.format(self.file['result'])

        # runcmd(self, 'run_fastq_dup', cmd)
        command = self.add_command("run_fastq_dup", cmd ,ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("succeed in running fastq_dup")
        elif command.return_code == 1:
            self.logger.info("return code: {}".format(command.return_code))
            self.add_state('memory_limit', 'memory is low!')
        elif command.return_code is None:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("succeed in running fastq_dup")
        else:
            self.set_error("fail to run fastq_dup")

    def set_output(self):
        self.option('result').set_path(self.file['result'])
