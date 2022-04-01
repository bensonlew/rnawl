# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest
import pandas as pd

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class KmSurvivalAgent(Agent):
    '''
    last_modify: 2020.06.29
    '''

    def __init__(self, parent):
        super(KmSurvivalAgent, self).__init__(parent)
        options = [
            {'name': 'km_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'risk_table', 'type': 'bool', 'default': True},
            {'name': 'conf', 'type': 'bool', 'default': True},
            {'name': 'km_pdf', 'type': 'outfile', 'format':'whole_transcriptome.common'}

        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(KmSurvivalAgent, self).end()


class KmSurvivalTool(Tool):
    def __init__(self, config):
        super(KmSurvivalTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.program = {
            'Rscript': 'program/R-3.3.1/bin/Rscript',
        }
        self.script = {
            'km_survival': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/km_survival.r')
        }
        self.file = {
            'km_pdf': os.path.join(self.output_dir, "km.pdf")
        }

    def run(self):
        super(KmSurvivalTool, self).run()
        self.run_km_survival()
        self.set_output()
        self.end()

    def run_km_survival(self):
        cmd = '{} {} '.format(self.program['Rscript'], self.script['km_survival'])
        cmd += ' -k {}'.format(self.option('km_table').path)
        cmd += ' -r {}'.format(self.option('risk_table'))
        cmd += ' -c {}'.format(self.option('conf'))
        cmd += ' -o {}'.format(self.file['km_pdf'])
        cmd_name = 'run_km_survival'
        command = self.add_command(cmd_name, cmd)
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
                self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704306")
        else:
            self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd), code = "33704307")



    def set_output(self):
        self.option('km_pdf').set_path(self.file['km_pdf'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'km_survival_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'tool_lab.km_survival',
            'instant': False,
            'options': {
                'km_table': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/km_survival/km_survival.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)


