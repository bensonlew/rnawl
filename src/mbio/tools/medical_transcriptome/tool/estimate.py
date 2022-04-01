# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import json
import os
import shutil
import unittest
import math
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class EstimateAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(EstimateAgent, self).__init__(parent)
        options = [
            {'name': 'exp_file', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'id', 'type': 'string', 'default': 'GeneSymbol'},
            {'name': 'platform', 'type': 'string', 'default': 'illumina'},
            {'name': 'species', 'type': 'string', 'default': 'Homo_sapiens'},
            {'name': 'estimate_score_gct', 'type': 'outfile', 'format': 'medical_transcriptome.common'},
            {'name': 'cluster_matrix', 'type': 'outfile', 'format': 'medical_transcriptome.common'}


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
        super(EstimateAgent, self).end()


class EstimateTool(Tool):
    def __init__(self, config):
        super(EstimateTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.program = {
            'python': 'program/Python/bin/python',
            'rscript': 'program/R-3.3.1/bin/Rscript',
        }
        self.script = {
            'estimate': os.path.join(self.config.PACKAGE_DIR, 'medical_transcriptome/estimate.py'),
        }
        self.file = {
            'exp_new': os.path.join(self.work_dir, 'exp.txt'),
            'estimate_score_gct': os.path.join(self.output_dir, 'estimate_score.gct'),
            'Homo_sapiens': os.path.join(self.config.SOFTWARE_DIR, 'database/gene_db/hsapiens.all_merged.xls'),
            'Mus_musculus': os.path.join(self.config.SOFTWARE_DIR, 'database/gene_db/mmusculus.all_merged.xls'),
            'Rattus_norvegicus': os.path.join(self.config.SOFTWARE_DIR, 'database/gene_db/rnorvegicus.all_merged.xls'),
            'cluster_matrix': os.path.join(self.output_dir, 'estimate_score_matrix.txt')
        }

    def run(self):
        super(EstimateTool, self).run()
        self.prepare()
        self.estimate()
        self.estimate_matrix()
        self.set_output()
        self.end()

    def prepare(self):
        with open(self.file['{}'.format(self.option('species'))], 'r') as f:
            gene_id_name = dict()
            for line in f.readlines():
                ensembl_id = line.strip().split('\t')[4]
                gene_name = line.strip().split('\t')[5]
                gene_id_name[ensembl_id] = gene_name
        exp = pd.read_table(self.option('exp_file').path, sep='\t', header=0, index_col=None)
        seq_id = exp['seq_id']
        seq_name = list()
        for i in seq_id:
            if i in gene_id_name:
                seq_name.append(gene_id_name[i])
            else:
                seq_name.append('-')
        exp['seq_id'] = seq_name
        exp = exp.drop(exp[exp['seq_id'] == '-'].index)
        exp = exp.drop_duplicates(['seq_id'])
        exp.to_csv(self.file['exp_new'], header=True, index=False, sep='\t')
    def estimate(self):
        cmd = '{} {}'.format(self.program['python'], self.script['estimate'])
        cmd += ' -e {}'.format(self.file['exp_new'])
        cmd += ' -i {}'.format(self.option('id'))
        cmd += ' -p {}'.format(self.option('platform'))
        cmd += ' -o {}'.format(self.output_dir)
        cmd_name = 'run_estimate'
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

    def estimate_matrix(self):
        es = pd.read_table(self.file['estimate_score_gct'], skiprows=2, sep='\t', header=0, index_col=0)
        es = es.drop(columns='Description')
        es_score = es.loc['ESTIMATEScore'].tolist()
        TumorPurity = [math.cos(0.6049872018+0.0001467884*x) for x in es_score]
        es.loc['TumorPurity'] = TumorPurity
        es.to_csv(self.file['cluster_matrix'], header=True, index=True, sep='\t')

    def set_output(self):
        self.option('estimate_score_gct').set_path(self.file['estimate_score_gct'])
        self.option('cluster_matrix').set_path(self.file['cluster_matrix'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'estimate{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'medical_transcriptome.tool.estimate',
            'instant': False,
            'options': {
                'exp_file': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/nomogram/test2medical/count_test.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
