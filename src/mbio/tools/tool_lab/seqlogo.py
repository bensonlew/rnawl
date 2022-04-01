# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'


import unittest
import os
import glob
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class SeqlogoAgent(Agent):
    def __init__(self, parent):
        super(SeqlogoAgent, self).__init__(parent)
        options = [
            {"name": "infile", "type": "infile", 'format': 'ref_rna_v2.common'},
            {"name": "matrix_type", "type": "string", "format": "none"},
            {'name': 'seq_type', 'type': 'string', 'default': 'nt'},
            {'name': 'unit', 'type': 'string', 'default': 'both'},
        ]
        self.add_option(options)
        self.step.add_steps("seqlogo")
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.seqlogo.start()
        self.step.update()

    def step_finish(self):
        self.step.seqlogo.finish()
        self.step.update()

    def check_options(self):
        if not self.option("infile").is_set:
            raise OptionError("必须设置输入序列/矩阵文件。")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        super(SeqlogoAgent, self).end()


class SeqlogoTool(Tool):
    def __init__(self, config):
        super(SeqlogoTool, self).__init__(config)
        self._version = "v1.0"
        self.program = {
            'python': os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/rna/miniconda3/bin/python3')
        }
        self.script = {
            'seqlogo': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/seqlogo.py')
        }

    def run(self):
        super(SeqlogoTool, self).run()
        self.run_seqlogo()
        self.set_output()
        self.end()

    def run_seqlogo(self):
        cmd = '{} {}'.format(self.program['python'], self.script['seqlogo'])
        cmd += ' -i {}'.format(self.option('infile').prop['path'])
        cmd += ' -m {}'.format(self.option('matrix_type'))
        cmd += ' -u {}'.format(self.option('unit'))
        cmd += ' -s {}'.format(self.option('seq_type'))
        cmd_name = 'draw_seqlogo'
        command = self.add_command(cmd_name, cmd, shell=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} 运行成功".format(cmd_name))
            else:
                self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd))
        else:
            self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd))

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        files = glob.glob(os.path.join(self.work_dir, '*.pdf'))
        for each in files:
            new_path = os.path.join(self.output_dir, os.path.basename(each))
            os.link(each, new_path)

