# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
#from mbio.packages.denovo_rna.assemble.trinity_stat import *
import os
import re
import shutil
import unittest
import pandas as pd


class BuscoBarAgent(Agent):
    """
    Trinity拼接
    author: 刘彬旭
    last_modify: 2017.11.9
    """
    def __init__(self, parent):
        super(BuscoBarAgent, self).__init__(parent)

        options = [
            {"name": "summary_result", "type": "infile", "format": "denovo_rna_v2.common"},
            {"name": "busco_result", "type": "outfile", "format": "denovo_rna_v2.common"},
            {'name': 'title', 'type': 'string'},
            {'name': 'xlab', 'type': 'string'},
        ]
        self.add_option(options)


    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        pass


    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(BuscoBarAgent, self).end()


class BuscoBarTool(Tool):
    def __init__(self, config):
        super(BuscoBarTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.program = {
            'python': 'program/Python/bin/python',
            'rscript': 'program/R-3.3.1/bin/Rscript',
        }
        self.script = {
            'busco_bar': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/busco_bar.r'),
        }
        self.file = {
            'busco_prepare': os.path.join(self.work_dir, 'busco_bar.txt'),
            'pdf_output': os.path.join(self.output_dir, 'busco.pdf')
        }

    def run(self):
        """
        运行
        :return:
        """
        super(BuscoBarTool, self).run()
        self.prepare()
        self.run_busco_bar()
        self.set_output()
        self.end()

    def prepare(self):
        busco = pd.read_table(self.option('summary_result').path, header=0, index_col=0, sep='\t')
        need_list = busco.index.tolist()[1:-1]
        busco_new = busco.loc[need_list]
        busco_new['spe'] = 'spe1'
        busco_new.to_csv(self.file['busco_prepare'], header=True, index=True, sep='\t')

    def run_busco_bar(self):
        cmd = '{} {} {} {} {} {}'.format(self.program['rscript'], self.script['busco_bar'], self.file['busco_prepare'],
                                   self.option('title'), self.option('xlab'), self.file['pdf_output'])
        cmd_name = 'run_busco_bar'
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
        pass

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'busco{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'tool_lab.busco',
            'instant': False,
            'options': {
                'fa': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/Busco/assemble_raw.fasta',
                'odb9': 'eukaryota_odb9',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)