# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest
import pandas as pd
import numpy as np

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class NineQuadrantAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(NineQuadrantAgent, self).__init__(parent)
        options = [
            {'name': 'trans_table', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'protein_table', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'id_trans_table', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'trans_number', 'type': 'int'},
            {'name': 'protein_number', 'type': 'int'},
            {'name': 'trans_fc', 'type': 'float', 'default': 2},
            {'name': 'protein_fc', 'type': 'float', 'default': 1.5},
            {'name': 'colour', 'type': 'string', 'default': 'Set1'},
            {'name': 'x_lab', 'type': 'string', 'default': 'log2_ratio_of_protein'},
            {'name': 'y_lab', 'type': 'string', 'default': 'log2_ratio_of_transcript'},
            {'name': 'log2fc_table', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'pdf_path', 'type': 'outfile', 'format': 'ref_rna_v2.common'}

        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '20G'

    def end(self):
        super(NineQuadrantAgent, self).end()


class NineQuadrantTool(Tool):
    def __init__(self, config):
        super(NineQuadrantTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.program = {
            'python': 'miniconda2/bin/python',
            'r': '/program/R-3.3.1/bin/Rscript'
        }
        self.script = {
            'nine_quadrant': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/nine_quadrant.r')
        }
        self.file = {
            'log2fc_table': os.path.join(self.output_dir, '{}'.format('log2fc.txt')),
            'pdf_path': os.path.join(self.output_dir, '{}'.format('nine_quadrant.pdf'))
        }

    def run(self):
        super(NineQuadrantTool, self).run()
        self.pre_nine_quadrant()
        self.run_nine_quadrant()
        self.set_output()
        self.end()

    def pre_nine_quadrant(self):
        if self.option('id_trans_table').is_set:
            id_trans_dict = dict()
            # with open(self.option('id_trans_table'), 'w') as t:
            with open(self.option('id_trans_table').prop['path'], 'r') as t:
                for i in t.readlines():
                    # protein_id, trans_id = i.split().strip('\t')
                    trans_id, protein_id= i.strip().split('\t')
                    id_trans_dict[protein_id] = trans_id

        df_trans = pd.read_table(self.option('trans_table').path, sep='\t', header=0, index_col=0)
        df_protein = pd.read_table(self.option('protein_table').path, sep='\t', header=0, index_col=0)
        df_log2 = pd.DataFrame(columns=['trans_log2fc', 'protein_log2fc'])
        trans_list = df_trans.columns.tolist()
        protein_list = df_protein.columns.tolist()
        trans_fc = trans_list[self.option('trans_number')-2]
        protein_fc = protein_list[self.option('protein_number')-2]


        if not self.option('id_trans_table').is_set:
            trans_log2fc = list()
            protein_log2fc = list()
            log2_index = list()
            for i in df_protein.index:
                if i in df_trans.index:
                    log2_index.append(i)
                    log2fc_p = df_protein.loc[i][protein_fc]
                    protein_log2fc.append(log2fc_p)
                    log2fc_t = df_trans.loc[i][trans_fc]
                    trans_log2fc.append(log2fc_t)
                else:
                    continue
            df_log2['trans_log2fc'] = trans_log2fc
            df_log2['protein_log2fc'] = protein_log2fc
            df_log2['id'] = log2_index
        else:
            trans_log2fc = list()
            protein_log2fc = list()
            protein_log2_index = list()
            trans_log2_index = list()
            for i in df_protein.index:
                # if id_trans_dict[i] in df_trans.index:
                if id_trans_dict.get(i) in df_trans.index:
                    protein_log2_index.append(i)
                    trans_log2_index.append(id_trans_dict[i])
                    log2fc_p = df_protein.loc[i][protein_fc]
                    protein_log2fc.append(log2fc_p)
                    log2fc_t = df_trans.loc[id_trans_dict[i]][trans_fc]
                    trans_log2fc.append(log2fc_t)
                else:
                    continue
            df_log2['trans_log2fc'] = trans_log2fc
            df_log2['protein_log2fc'] = protein_log2fc
            df_log2['protein_id'] = protein_log2_index
            df_log2['trans_id'] = trans_log2_index

        group = list()
        for i in df_log2.index:
            if df_log2.loc[i]['trans_log2fc'] >= np.log2(self.option('trans_fc')) and df_log2.loc[i]['protein_log2fc'] >= np.log2(self.option('protein_fc')):
                group.append('quadrant3')
            if df_log2.loc[i]['trans_log2fc'] >= np.log2(self.option('trans_fc')) and -np.log2(self.option('protein_fc')) < df_log2.loc[i]['protein_log2fc'] < np.log2(self.option('protein_fc')):
                group.append('quadrant2')
            if df_log2.loc[i]['trans_log2fc'] >= np.log2(self.option('trans_fc')) and df_log2.loc[i]['protein_log2fc'] <= -np.log2(self.option('protein_fc')):
                group.append('quadrant1')
            if df_log2.loc[i]['trans_log2fc'] <= -np.log2(self.option('trans_fc')) and df_log2.loc[i]['protein_log2fc'] >= np.log2(self.option('protein_fc')):
                group.append('quadrant9')
            if df_log2.loc[i]['trans_log2fc'] <= -np.log2(self.option('trans_fc')) and -np.log2(self.option('protein_fc')) < df_log2.loc[i]['protein_log2fc'] < np.log2(self.option('protein_fc')):
                group.append('quadrant8')
            if df_log2.loc[i]['trans_log2fc'] <= -np.log2(self.option('trans_fc')) and df_log2.loc[i]['protein_log2fc'] <= -np.log2(self.option('protein_fc')):
                group.append('quadrant7')
            if -np.log2(self.option('trans_fc')) < df_log2.loc[i]['trans_log2fc'] < np.log2(self.option('trans_fc')) and df_log2.loc[i]['protein_log2fc'] >= np.log2(self.option('protein_fc')):
                group.append('quadrant6')
            if -np.log2(self.option('trans_fc')) < df_log2.loc[i]['trans_log2fc'] < np.log2(self.option('trans_fc')) and -np.log2(self.option('protein_fc')) < df_log2.loc[i]['protein_log2fc'] < np.log2(self.option('protein_fc')):
                group.append('quadrant5')
            if -np.log2(self.option('trans_fc')) < df_log2.loc[i]['trans_log2fc'] < np.log2(self.option('trans_fc')) and df_log2.loc[i]['protein_log2fc'] <= -np.log2(self.option('protein_fc')):
                group.append('quadrant4')
        df_log2['group'] = group
        if not self.option('id_trans_table').is_set:
            df_log2_new = df_log2[['id', 'trans_log2fc', 'protein_log2fc', 'group']]
        else:
            df_log2_new = df_log2[['protein_id', 'trans_id', 'trans_log2fc', 'protein_log2fc', 'group']]
        df_log2_new.to_csv(self.file['log2fc_table'], sep='\t', header=True, index=False)

    def run_nine_quadrant(self):
        cmd = '{} {} {} {} {} {} {} {} {}'.format(self.program['r'], self.script['nine_quadrant'], self.file['log2fc_table'],
                                                  self.option('protein_fc'), self.option('trans_fc'), self.file['pdf_path'],
                                                  self.option('colour'), self.option('x_lab'), self.option('y_lab'))
        cmd_name = 'run_nine_quadrant'
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
        self.option('log2fc_table').set_path(self.file['log2fc_table'])
        self.option('pdf_path').set_path(self.file['pdf_path'])

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'nine_quadrant_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'tool_lab.nine_quadrant',
            'instant': False,
            'options': {
                'trans_table': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/nine/to-vs-t1.genes.annot.txt',
                'protein_table': '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/tool_lab/nine/po-vs-p1.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)


