# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import os,glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
from mbio.packages.whole_transcriptome.utils import runcmd
import unittest
import pandas as pd
import numpy as np
from sklearn import decomposition, preprocessing


class RandomreadsAgent(Agent):
    def __init__(self, parent):
        super(RandomreadsAgent, self).__init__(parent)
        options = [
            {'name': 'ref_fa', 'type': 'string'},
            {'name': 'paired', 'type': 'string'},
            {'name': 'split', 'type': 'string'},
            {'name': 'reads', 'type': 'int'},
            {'name': 'reads_length', 'type': 'int'},
            {'name': 'pacbio', 'type': 'string'},
            {'name': 'mininsert', 'type': 'int'},
            {'name': 'maxinsert', 'type': 'int'},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = "60G"

    def end(self):
        super(RandomreadsAgent, self).end()


class RandomreadsTool(Tool):
    def __init__(self, config):
        super(RandomreadsTool, self).__init__(config)
        java_dir = os.path.join(self.config.SOFTWARE_DIR, 'program/sun_jdk1.8.0/bin')
        self.set_environ(PATH=java_dir)
        self.program = {
            'python': 'miniconda2/bin/python',
        }
        self.script = {
            'randomreads': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/Randomreads.py'),
            'bbmap_randomreads': os.path.join(self.config.SOFTWARE_DIR, "bioinfo/tool_lab/bbmap/randomreads.sh"),
            'bbmap_reformat': os.path.join(self.config.SOFTWARE_DIR, "bioinfo/tool_lab/bbmap/reformat.sh")
        }
        self.file = {
            'output_file': os.path.join(self.output_dir, os.path.basename(self.option('ref_fa')).split('.')[0]) + '.fq.gz',
            'out1': os.path.join(self.output_dir, os.path.basename(self.option('ref_fa')).split('.')[0]) + '.R1.fq.gz',
            'out2': os.path.join(self.output_dir, os.path.basename(self.option('ref_fa')).split('.')[0]) + '.R2.fq.gz',
        }
    def run(self):
        """
        运行
        :return:
        """
        super(RandomreadsTool, self).run()
        self.run_randomreads()
        if self.option('split') == 't':
            self.run_randomreads_split()
        self.set_output()
        self.end()

    def run_randomreads_split(self):
        cmd = '{}'.format(self.script['bbmap_reformat'])
        cmd += ' in={}'.format(self.file['output_file'])
        cmd += ' out1={}'.format(self.file['out1'])
        cmd += ' out2={}'.format(self.file['out2'])
        cmd += ' overwrite=t'
        cmd += ' addcolon'
        cmd_name = 'run_randomreads_split'
        command = self.add_command(cmd_name, cmd, shell=True)
        command.run()
        self.wait()
        if command.return_code in [0]:
            self.logger.info("{} Finished successfully".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} Finished successfully".format(cmd_name))
            else:
                self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="33704402")
        else:
            self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="33704403")
    def run_randomreads(self):
        cmd = '{}'.format(self.script['bbmap_randomreads'])
        if self.option('split') == 'f':
            cmd += ' ref={}'.format(self.option('ref_fa'))
            cmd += ' out={}'.format(self.file['output_file'])
            cmd += ' length={}'.format(self.option('reads_length'))
            cmd += ' paired={}'.format(self.option('paired'))
            cmd += ' mininsert={}'.format(self.option('mininsert'))
            cmd += ' maxinsert={}'.format(self.option('maxinsert'))
            cmd += ' reads={}'.format(self.option('reads'))
            cmd += ' pacbio={}'.format(self.option('pacbio'))
            cmd += ' overwrite=t'
        if self.option('split') == 't':
            cmd += ' ref={}'.format(self.option('ref_fa'))
            cmd += ' out={}'.format(self.file['output_file'])
            cmd += ' length={}'.format(self.option('reads_length'))
            cmd += ' paired={}'.format(self.option('paired'))
            cmd += ' mininsert={}'.format(self.option('mininsert'))
            cmd += ' maxinsert={}'.format(self.option('maxinsert'))
            cmd += ' reads={}'.format(self.option('reads'))
            cmd += ' pacbio={}'.format(self.option('pacbio'))
        cmd_name = 'run_randomreads'
        command = self.add_command(cmd_name, cmd, shell=True, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code in [0]:
            self.logger.info("{} Finished successfully".format(cmd_name))
        elif command.return_code is None or command.return_code == 127:
            self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} Finished successfully".format(cmd_name))
            else:
                self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="33704402")
        else:
            self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="33704403")

    # def run_randomreads(self):
    #     cmd = '{} {}'.format(self.program['python'], self.script['randomreads'])
    #     cmd += ' -b {}'.format(self.script['bbmap'])
    #     cmd += ' -r {}'.format(self.option('ref_fa'))
    #     cmd += ' -o {}'.format(self.file['output_file'])
    #     cmd += ' -l {}'.format(self.option('reads_length'))
    #     cmd += ' -pair {}'.format(self.option('paired'))
    #     cmd += ' -mininsert {}'.format(self.option('mininsert'))
    #     cmd += ' -maxinsert {}'.format(self.option('maxinsert'))
    #     cmd += ' -reads {}'.format(self.option('reads'))
    #     cmd += ' -pacbio {}'.format(self.option('pacbio'))
    #     cmd += ' -split {}'.format(self.option('split'))
    #     cmd_name = 'run_randomreads'
    #     command = self.add_command(cmd_name, cmd)
    #     command.run()
    #     self.wait()
    #     if command.return_code in [0]:
    #         self.logger.info("{} Finished successfully".format(cmd_name))
    #     elif command.return_code is None:
    #         self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_name))
    #         command.rerun()
    #         self.wait()
    #         if command.return_code is 0:
    #             self.logger.info("{} Finished successfully".format(cmd_name))
    #         else:
    #             self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="33704402")
    #     else:
    #         self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="33704403")



    def set_output(self):
        # self.option('translation_file').set_path(self.file['output_file'])
        pass

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "heatmap_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.heatmap",
            "options": dict(
                otutable="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/cluster/otu_table.xls",
                log_change="None",
                scale='zscore'

            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
