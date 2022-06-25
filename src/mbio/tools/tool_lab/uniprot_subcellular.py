# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import json
import os
import shutil
import unittest
import datetime
import time
import pickle
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class UniprotSubcellularAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(UniprotSubcellularAgent, self).__init__(parent)
        options = [
            {'name': 'uniprot_id', 'type': 'string'},
            {'name': 'method', 'type': 'string', 'default': 'uniprot_api'},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'
    def end(self):
        super(UniprotSubcellularAgent, self).end()


class UniprotSubcellularTool(Tool):
    def __init__(self, config):
        super(UniprotSubcellularTool, self).__init__(config)
        self.program = {
            'python': '/miniconda2/bin/python',
        }
        self.script = {
            'uniprot_api': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/Uniprot.py'),
        }
        self.file = {
            'subcellular_location': os.path.join(self.output_dir, 'subcellular_location.xlsx'),
        }

    def run(self):
        super(UniprotSubcellularTool, self).run()
        self.uniprot_id = ';'.join(self.option('uniprot_id').strip().split('\n'))
        time_file = self.config.PACKAGE_DIR + '/tool_lab/uniprot.time'
        if not os.path.exists(time_file):
            with open(time_file, 'w') as tw:
                pickle.dump(datetime.datetime.now(), tw)
        wait_start = datetime.datetime.now()
        while 1:
            tr = open(time_file)
            last_time = pickle.load(tr)
            now_time = datetime.datetime.now()
            if not last_time:
                if (now_time - wait_start).seconds < 720:
                    time.sleep(60)
                    continue
                else:
                    break

            # if (now_time-last_time).seconds < 720:
            if (now_time-last_time).seconds < 20:
                self.logger.info('需要睡1分钟')
                time.sleep(60)
            else:
                break
        with open(time_file, 'w') as tw:
            pickle.dump(datetime.datetime.now(), tw)
            # pickle.dump('', tw)
        self.run_uniprot_api()
        self.set_output()
        self.end()

    def run_uniprot_api(self):
        cmd = '{} {}'.format(self.program['python'], self.script['uniprot_api'])
        cmd += ' -uniprot_id "{}"'.format(self.uniprot_id)
        cmd += ' -output {}'.format(self.file['subcellular_location'])
        cmd += ' -method {}'.format(self.option('method'))
        cmd_name = 'run_uniprot'
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
                self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="33704402")
        else:
            self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="33704403")



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
            'id': 'protein_picture_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'tool_lab.protein_picture',
            'instant': False,
            'options': {
                'picture': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/protein_picture/Dl09YYJ_bjhb1_1580_Result.PNG',
                'if_black': 'no',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
