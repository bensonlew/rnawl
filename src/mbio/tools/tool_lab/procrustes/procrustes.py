# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import glob
import unittest
import re


class ProcrustesAgent(Agent):
    """
    PCoA
    """
    def __init__(self, parent):
        super(ProcrustesAgent, self).__init__(parent)
        options = [
            {'name': 'coord_ref', 'type': 'infile', 'default': 'None', 'format': 'metabolome.coord_matrices',
             'required': True},  # 只能有一个
            {'name': 'coord_query', 'type': 'infile', 'default': 'None', 'format': 'metabolome.coord_matrices',
             'required': True},  # 如有多个可逗号隔开，每个均会与ref进行计算
            {'name': 'random', 'type': 'int', 'default': '1000', 'format': 'None'},
            {'name': 'out', 'type': 'string', 'default': self.work_dir, 'format': 'None'}
        ]
        self.add_option(options)
        self.step.add_steps('procrustes')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.procrustes.start()
        self.step.update()

    def step_end(self):
        self.step.procrustes.finish()
        self.step.update()

    def check_options(self):
        if not self.option('coord_ref').is_set:
            raise OptionError('必须提供输入ref矩阵表')
        if not self.option('coord_query').is_set:
            raise OptionError('必须提供输入query矩阵表')
        return True

    def set_resource(self):
        """
        运行所需资源
        """
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        super(ProcrustesAgent, self).end()


class ProcrustesTool(Tool):
    def __init__(self, config):
        super(ProcrustesTool, self).__init__(config)
        self._version = "1.0"
        self.gcc = os.path.join(self.config.SOFTWARE_DIR, 'gcc/5.1.0/bin')
        self.gcc_lib = os.path.join(self.config.SOFTWARE_DIR, 'gcc/5.1.0/lib64')
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.program = {
            'python': "miniconda2/bin/python",
        }
        self.script = {
            'procrustes': os.path.join(self.config.SOFTWARE_DIR, "miniconda2/bin/transform_coordinate_matrices.py"),
        }

    def run_procrustes(self):
        cmd = '{} {} '.format(self.program['python'], self.script['procrustes'])
        cmd += '-i {},{} '.format(self.option("coord_ref").prop["path"], self.option("coord_query").prop["path"])
        cmd += '-o {} '.format(self.option("out"))
        cmd += '-r {} '.format(self.option("random"))
        cmd_name = 'procrustes'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} Finished successfully".format(cmd_name))
            else:
                self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))
        else:
            self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))

    def set_output(self):
        target_files = glob.glob(self.option("out") + '/*transformed_reference.txt')
        target_files += glob.glob(self.option("out") + '/*transformed_q*.txt')
        target_files += glob.glob(self.option("out") + '/procrustes_results.txt')
        for each in target_files:
            name = os.path.basename(each)
            link = os.path.join(self.output_dir, name)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(ProcrustesTool, self).run()
        self.run_procrustes()
        self.set_output()
        self.end()

