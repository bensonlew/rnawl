# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import shutil
import unittest

import numpy as np
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool

from mbio.packages.whole_transcriptome.utils import runcmd


class ExpMakeAgent(Agent):
    '''
    last_modify: 2019.10.22
    '''

    def __init__(self, parent):
        super(ExpMakeAgent, self).__init__(parent)
        options = [
            {'name': 'exp_class', 'type': 'string', 'default': 'mlc'},
            {'name': 't_tpm', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'g_tpm', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 't_fpkm', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'g_fpkm', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 't_count', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'g_count', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'annot_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'relat_table', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'c_rpm', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'c_detail', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'c_count', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 's_count', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'count_dir', 'type': 'outfile', 'format': 'whole_transcriptome.common_dir'},
        ]
        self.add_option(options)

    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '20G'

    def end(self):
        super(ExpMakeAgent, self).end()


class ExpMakeTool(Tool):
    def __init__(self, config):
        super(ExpMakeTool, self).__init__(config)
        self.program = {
            'python': 'program/Python/bin/python'
        }
        self.script = {
            'exp_annot': os.path.join(self.config.PACKAGE_DIR, 'whole_transcriptome/formation/exp_annot.py')
        }
        self.dir = {
            'tpm': os.path.join(self.work_dir, 'tpm'),
            'fpkm': os.path.join(self.work_dir, 'fpkm'),
            'count': os.path.join(self.work_dir, 'count'),
            'm_output': os.path.join(self.output_dir, 'mrna'),
            'l_output': os.path.join(self.output_dir, 'lncrna'),
            'c_output': os.path.join(self.output_dir, 'circrna'),
            's_output': os.path.join(self.output_dir, 'smallrna'),
        }

    def run(self):
        super(ExpMakeTool, self).run()
        if 'm' in self.option('exp_class'):
            self.run_exp_annot_tpm()
            self.run_exp_annot_fpkm()
        self.collect_count_matrix()
        self.set_output()
        self.end()

    def run_exp_annot_tpm(self):
        if os.path.isdir(self.dir['tpm']):
            shutil.rmtree(self.dir['tpm'])
        os.makedirs(self.dir['tpm'])
        cmd = '{} {}'.format(self.program['python'], self.script['exp_annot'])
        cmd += ' -t {}'.format(self.option('t_tpm').path)
        cmd += ' -g {}'.format(self.option('g_tpm').path)
        cmd += ' -a {}'.format(self.option('annot_table').path)
        cmd += ' -r {}'.format(self.option('relat_table').path)
        cmd += ' -o {}'.format(self.dir['tpm'])
        runcmd(self, 'run_exp_annot_tpm', cmd, block=False)

    def run_exp_annot_fpkm(self):
        if os.path.isdir(self.dir['fpkm']):
            shutil.rmtree(self.dir['fpkm'])
        os.makedirs(self.dir['fpkm'])
        cmd = '{} {}'.format(self.program['python'], self.script['exp_annot'])
        cmd += ' -t {}'.format(self.option('t_fpkm').path)
        cmd += ' -g {}'.format(self.option('g_fpkm').path)
        cmd += ' -a {}'.format(self.option('annot_table').path)
        cmd += ' -r {}'.format(self.option('relat_table').path)
        cmd += ' -o {}'.format(self.dir['fpkm'])
        runcmd(self, 'run_exp_annot_fpkm', cmd)

    def collect_count_matrix(self):
        if os.path.isdir(self.dir['count']):
            shutil.rmtree(self.dir['count'])
        os.makedirs(self.dir['count'])
        file_num = 0
        for level in ['T', 'G']:
            if 'm' in self.option('exp_class'):
                count_df = pd.read_table(self.option('{}_count'.format(level.lower())).path)
                if 'rna_type' in count_df.columns and 'is_new' in count_df.columns:
                    count_df = count_df.iloc[:, :-2]
                for idx in count_df.columns[1:]:
                    count_df[idx] = count_df[idx].apply(round).astype(np.int32)
                else:
                    count_df.to_csv(os.path.join(self.dir['count'], '{}.reads.txt'.format(level)), sep='\t',
                                    index=False)
                file_num += 1
        else:
            if 'c' in self.option('exp_class'):
                shutil.copy(self.option('c_count').path, os.path.join(self.dir['count'], 'C.reads.txt'))
                file_num += 1
            if 's' in self.option('exp_class'):
                shutil.copy(self.option('s_count').path, os.path.join(self.dir['count'], 'S.reads.txt'))
                file_num += 1
            self.logger.info('succeed in exporting {} count files in {}'.format(file_num, self.dir['count']))

    def set_output(self):
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        if 'm' in self.option('exp_class'):
            os.mkdir(self.dir['m_output'])
            os.link(os.path.join(self.dir['tpm'], 'T.mRNA.txt'), os.path.join(self.dir['m_output'], 'T.tpm.txt'))
            os.link(os.path.join(self.dir['tpm'], 'G.mRNA.txt'), os.path.join(self.dir['m_output'], 'G.tpm.txt'))
            os.link(os.path.join(self.dir['fpkm'], 'T.mRNA.txt'), os.path.join(self.dir['m_output'], 'T.fpkm.txt'))
            os.link(os.path.join(self.dir['fpkm'], 'G.mRNA.txt'), os.path.join(self.dir['m_output'], 'G.fpkm.txt'))
        if 'l' in self.option('exp_class'):
            os.mkdir(self.dir['l_output'])
            os.link(os.path.join(self.dir['tpm'], 'T.lncRNA.txt'), os.path.join(self.dir['l_output'], 'T.tpm.txt'))
            os.link(os.path.join(self.dir['tpm'], 'G.lncRNA.txt'), os.path.join(self.dir['l_output'], 'G.tpm.txt'))
            os.link(os.path.join(self.dir['fpkm'], 'T.lncRNA.txt'), os.path.join(self.dir['l_output'], 'T.fpkm.txt'))
            os.link(os.path.join(self.dir['fpkm'], 'G.lncRNA.txt'), os.path.join(self.dir['l_output'], 'G.fpkm.txt'))
        self.prepare_longrna_data()
        if 'c' in self.option('exp_class'):
            os.mkdir(self.dir['c_output'])
            c_exp_table = os.path.join(self.dir['c_output'], 'T.rpm.txt')
            shutil.copy(self.option('c_rpm').path, c_exp_table)
            c_exp_df = pd.read_table(c_exp_table, index_col=0)
            c_exp_df.index.name = 'transcript_id'
            c_det_df = pd.read_table(self.option('c_detail').path, index_col=0)
            if 'circbase' in c_det_df.columns:
                c_base_dict = c_det_df['circbase'].to_dict()
            else:
                c_base_dict = dict()
            c_gid_df = pd.read_table(self.option('c_detail').path, header=0, names=['transcript_id', 'gene_id'],
                                     index_col=0, usecols=[0, 1])
            c_exp_df = c_exp_df.join(c_gid_df)
            kinds = list()
            for transcript_id in c_exp_df.index:
                if c_base_dict.get(transcript_id) == 'yes':
                    kinds.append('ref')
                else:
                    kinds.append('new')
            c_exp_df['level'] = 'T'
            c_exp_df['category'] = 'circRNA'
            c_exp_df['kind'] = kinds
            c_exp_df.reset_index(inplace=True)
            c_exp_df = c_exp_df.drop_duplicates()
            c_exp_df.to_csv(c_exp_table, sep='\t', index=False)
        if 's' in self.option('exp_class'):
            pass
        shutil.copytree(self.dir['count'], os.path.join(self.output_dir, 'count'))
        self.option('count_dir').set_path(os.path.join(self.output_dir, 'count'))

    def prepare_longrna_data(self):
        longrna_dir = os.path.join(self.output_dir, 'longrna')
        if os.path.isdir(longrna_dir):
            shutil.rmtree(longrna_dir)
        if 'l' in self.option('exp_class'):
            os.makedirs(longrna_dir)
            for way, level in (('tpm', 'T'), ('tpm', 'G'), ('fpkm', 'T'), ('fpkm', 'G')):
                output_path = os.path.join(longrna_dir, '{}.{}.txt'.format(level, way))
                m_df = pd.read_table(os.path.join(self.dir[way], '{}.mRNA.txt'.format(level)))
                l_df = pd.read_table(os.path.join(self.dir[way], '{}.lncRNA.txt'.format(level)))
                columns = l_df.columns.drop('gene_id') if level == 'T' else l_df.columns
                df = pd.concat([m_df.reindex(columns, axis=1), l_df.reindex(columns, axis=1)])
                df.rename({'{}'.format({'T': 'transcript_id', 'G': 'gene_id'}[level]): 'seq_id'}, axis=1)
                df.to_csv(output_path, sep='\t', index=False)
                self.logger.debug('succeed in exporting {}'.format(output_path))
        else:
            shutil.copytree(self.dir['m_output'], longrna_dir)
        self.logger.info('succeed in preparing longRNA expression matrix')


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'exp_make_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.formation.exp_make',
            'instant': False,
            'options': {
                'exp_class': 'mlc',
                't_tpm': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/LargeGush/output/filter_by_express/classifyquant/T.tpm.txt',
                'g_tpm': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/LargeGush/output/filter_by_express/classifyquant/G.tpm.txt',
                't_fpkm': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/LargeGush/output/filter_by_express/classifyquant/T.fpkm.txt',
                'g_fpkm': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/LargeGush/output/filter_by_express/classifyquant/G.fpkm.txt',
                't_count': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/LargeGush/output/filter_by_express/classifyquant/T.count.txt',
                'g_count': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/LargeGush/output/filter_by_express/classifyquant/G.count.txt',
                'annot_table': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/Annotation/output/allannot_class/all_annot.xls',
                'relat_table': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/LargeGush/output/filter_by_express/filtered_file/trans_type.xls',
                'c_rpm': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/CircBrush/output/RPM.txt',
                'c_detail': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/CircBrush/output/detail.txt',
                'c_count': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/CircBrush/output/count.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
