# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest
import datetime
import time
import pickle

from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class DiffStringPicturesAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(DiffStringPicturesAgent, self).__init__(parent)
        options = [
            {'name': 'string_xml', 'type': 'infile', 'format': 'itraq_and_tmt.common'},
            {'name': 'gene_list', 'type': 'string'},
            {'name': 'identity', 'type': 'float', 'default': 98},
            {'name': 'max_num', 'type': 'int', 'default': 300},
            {'name': 'species', 'type': 'int', 'default': 0},
            {'name': 'useblast', 'type': 'string', 'default': 'no'},

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
        super(DiffStringPicturesAgent, self).end()


class DiffStringPicturesTool(Tool):
    def __init__(self, config):
        super(DiffStringPicturesTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.program = {
            'python': 'miniconda2/bin/python',
        }
        self.script = {
            'string': self.config.PACKAGE_DIR + '/tool_lab/get_string_picture.py',
        }
        self.file = {
            'cox_table': os.path.join(self.output_dir, 'cox.txt'),
            'cox_graph': os.path.join(self.output_dir, 'cox.png')
        }

    def run(self):
        super(DiffStringPicturesTool, self).run()
        time_file = self.config.PACKAGE_DIR + '/tool_lab/string.time'
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
                    time.sleep(300)
                    continue
                else:
                    break

            # if (now_time-last_time).seconds < 720:
            if (now_time-last_time).seconds < 20:
                self.logger.info('需要睡5分钟')
                time.sleep(300)
            else:
                break
        with open(time_file, 'w') as tw:
            pickle.dump(datetime.datetime.now(), tw)
            # pickle.dump('', tw)
        self.pre_string()
        self.string_pictures()
        self.set_output()
        self.end()

    def pre_string(self):
        gene_list = self.option('gene_list').strip('\n').split('\n')
        with open(os.path.join(self.work_dir, 'protein.list'), 'w') as p:
            for i in gene_list:
                p.write(i + '\n')

    def string_pictures(self):

        cmd = '{} {}'.format(self.program['python'], self.script['string'])
        cmd += ' -specie ' + str(self.option('species'))
        cmd += ' -list_path ' + self.work_dir
        if self.option('useblast') == 'yes':
            cmd += ' -vsstring ' + self.option('string_xml').prop['path']
        cmd += ' -identity ' + str(self.option('identity'))
        cmd += ' -max_num ' + str(self.option('max_num'))
        cmd += ' -useblast ' + self.option('useblast')
        cmd += ' -out ' + self.output_dir
        cmd_name = 'run_string'
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
                self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "33704402")
        else:
            self.set_error("%s Failed. >>>%s", variables = (cmd_name, cmd), code = "33704403")

    def set_output(self):
        # self.option('cox_table').set_path(self.file['cox_table'])
        # self.option('cox_graph').set_path(self.file['cox_graph'])
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
            'id': 'STING{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'tool_lab.diff_string_pictures',
            'instant': False,
            'options': {
                'gene_list': 'trpA',
                'species': '9606',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
