# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest
import pandas as pd

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError
from biocluster.config import Config


class IdnrAgent(Agent):
    '''
    last_modify: 2020.07.23
    '''

    def __init__(self, parent):
        super(IdnrAgent, self).__init__(parent)
        options = [
            {'name': 'species', 'type': 'string'},
            {'name': 'id_file', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'id2nr', 'type': 'outfile', 'format': 'whole_transcriptome.common'}

        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(IdnrAgent, self).end()


class IdnrTool(Tool):
    def __init__(self, config):
        super(IdnrTool, self).__init__(config)
        self._project_type = "ref_rna_v2"
        self.client = Config().get_mongo_client(mtype=self._project_type)
        self.db = self.client[Config().get_mongo_dbname(self._project_type)]
        self.db1 = Config().get_mongo_client(mtype=_project_type, dydb_forbid=True)[
            Config().get_mongo_dbname(_project_type, dydb_forbid=True)]

        self.connect = self.db1['sg_genome_db']
        self.record = self.connect.find_one({'name': self.option('species')})
        self.anno = self.record['anno_path_v2']

        self.program = {
            'python': 'miniconda2/bin/python'
        }
        self.script = {
            'id2nr': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/id2nr.py')
        }
        self.file = {
            'id2nr': os.path.join(self.output_dir, "description.txt"),
            'blast_nr': os.path.join(Config().SOFTWARE_DIR, 'database/Genome_DB_finish', self.anno, 'annot_mapdb/nr', 'blast.xml')
        }

    def run(self):
        super(IdnrTool, self).run()
        self.run_id2nr()
        self.set_output()
        self.end()

    def run_id2nr(self):

        cmd = '{} {} -xml {} -id {} -o {}'.format(self.program['python'], self.script['id2nr'], self.file['blast_nr'],
                                                  self.option('id_file').path, self.file['id2nr'])
        cmd_name = 'run_id2nr'
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
        self.option('id2nr').set_path(self.file['id2nr'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'id2nr_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'tool_lab.idnr',
            'instant': False,
            'options': {
                'species': 'Mus_musculus',
                'id_file': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/id2nr/test.txt'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)


