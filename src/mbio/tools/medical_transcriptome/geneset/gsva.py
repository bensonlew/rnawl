# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import os
import unittest
import pandas as pd
import json
from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.whole_transcriptome.utils import runcmd
from biocluster.core.exceptions import OptionError


class GsvaAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(GsvaAgent, self).__init__(parent)
        options = [
            {'name': 'matrix', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {'name': 'gmx', 'type': 'infile', 'format': 'ref_rna_v2.geneset_gmt'},
            {'name': 'geneset_source', 'type': 'string', 'default': 'msigdb'},
            {'name': 'genes_detail', 'type': 'infile', 'format': 'ref_rna_v2.common'},
            {"name": 'method', 'type': 'string', 'default': 'gsva'},
            {'name': 'kcdf', 'type': 'string', 'default': 'Gaussian'},
            {'name': 'min_num', 'type': 'float', 'default': 10},
            {'name': 'max_num', 'type': 'float', 'default': 500},
            {'name': 'es_method', 'type': 'string', 'default': 'diff'},
            {'name': 'es_exp_table', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'gsva', 'type': 'outfile', 'format': 'ref_rna_v2.common'}
        ]
        self.add_option(options)

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        for name in ('matrix', 'gmx'):
            if not self.option(name).is_set:
                raise OptionError("缺少 %s 输入文件" , variables=( name), code="33710602")
    def set_resource(self):
        self._cpu = 1
        self._memory = '20G'

    def end(self):
        super(GsvaAgent, self).end()


class GsvaTool(Tool):
    def __init__(self, config):
        super(GsvaTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python',
            'Rscript': 'program/R-3.3.1/bin/Rscript'
        }
        self.script = {
            'gsva': os.path.join(self.config.PACKAGE_DIR, 'medical_transcriptome/gsva.r')
        }
        self.file = {
            'json': os.path.join(self.work_dir, 'input.json'),
            'ensembl2symbol': os.path.join(self.work_dir, 'ensembl2symbol.txt'),
            'matrix2gsea': os.path.join(self.work_dir, '{}'.format('GSVA.txt')),
            'biomart': os.path.join(self.config.SOFTWARE_DIR,
                                    'database/Genome_DB_finish/vertebrates/Mus_musculus/GRCm38_Ensembl_96/biomart/biomart.txt'),
            'es_exp_table': os.path.join(self.output_dir, 'gsva_table.txt')
        }
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')

    def run(self):
        super(GsvaTool, self).run()
        self.matrix2gsea()
        self.run_pregsva()
        self.run_gsva()

        self.set_output()
        self.end()

    def matrix2gsea(self):
        matrix = pd.read_table(self.option('matrix').path,header=0,index_col=0,sep='\t')
        symbol = pd.read_table(self.option('genes_detail').path,header=0,index_col=None,sep='\t')
        # symbol = pd.read_table(self.file['ensembl2symbol'],header=0,index_col=None,sep='\t')
        # symbol.dropna(subset=['gene_name'],inplace=True)
        if self.option('geneset_source').lower() == 'msigdb':
            ensembl_id = symbol['gene_id']
            id_list = ensembl_id.tolist()
            matrix_1 = matrix.loc[id_list]
            ensemble_symbol_dict = dict()
            with open(self.option('genes_detail').path, 'r') as e2s:
                for line in e2s.readlines():
                    if 'gene_id' == line.strip().split('\t')[0]:
                        continue
                    gene_id, gene_name, description =line.strip().split('\t')
                    ensemble_symbol_dict[gene_id] = [gene_name, description]
            a = list()
            # c = list()
            for i in matrix_1.index.tolist():
                a.append(ensemble_symbol_dict[i][0])
                # c.append(ensemble_symbol_dict[i][1])
            b = [i.upper() for i in a]
            matrix_1['symble'] = b
            # matrix_1['description'] = c
            matrix_2 = matrix_1.groupby('symble')
            matrix_3 = matrix_2.max()
            matrix_3.index.name = 'geneid'
            col_list = matrix_3.columns.tolist()
            # col_list.remove('description')
            # col_list.insert(0,'description')
            # matrix_3['description'] = [ensemble_symbol_dict[l][1] for l in matrix_3.index.tolist()]
            matrix_4 = matrix_3[col_list]
            matrix_4.to_csv(self.file['matrix2gsea'], index=True, sep='\t')
        else:
            matrix.index.name = 'geneid'
            matrix.to_csv(self.file['matrix2gsea'], index=True, sep='\t')

    def run_pregsva(self):
        if self.option('es_method') == 'diff':
            es_diff = True
        if self.option('es_method') == 'max':
            es_diff = False

        json.dump({
            'matrix': self.file['matrix2gsea'],
            'gmx': self.option('gmx').path,
            'method': self.option('method'),
            'kcdf': self.option('kcdf'),
            'min_num': self.option('min_num'),
            'max_num': self.option('max_num'),
            'es_method': es_diff,
            'output': self.file['es_exp_table']
        }, open(self.file['json'], 'w'), indent=4)

    def run_gsva(self):
        cmd = '{} {} -i {}'.format(self.program['Rscript'], self.script['gsva'], self.file['json'])
        cmd_name = 'run_gsva'
        runcmd(self, cmd_name, cmd)

    def set_output(self):
        self.option('es_exp_table').set_path(self.file['es_exp_table'])
        self.option('gsva').set_path(self.file['matrix2gsea'])
class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'gsva_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'tool_lab.gsva',
            'instant': False,
            'options': {
                'matrix':'/mnt/ilustre/users/sanger-dev/workspace/20200612/GenesetGsea_tsg_37656_4600_807/gsea.txt',
                "gmx" : '/mnt/ilustre/users/sanger-dev/workspace/20200612/GenesetGsea_tsg_37656_4600_807/gsea.gmt',

            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)


