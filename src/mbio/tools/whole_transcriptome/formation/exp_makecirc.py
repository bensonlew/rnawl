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


class ExpMakecircAgent(Agent):
    '''
    last_modify: 2019.10.22
    '''

    def __init__(self, parent):
        super(ExpMakecircAgent, self).__init__(parent)
        options = [

            {'name': 'c_rpm', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'c_detail', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'c_count', 'type': 'infile', 'format': 'whole_transcriptome.common'},
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
        super(ExpMakecircAgent, self).end()


class ExpMakecircTool(Tool):
    def __init__(self, config):
        super(ExpMakecircTool, self).__init__(config)

        self.dir = {

            'count': os.path.join(self.work_dir, 'count'),


            'c_output': os.path.join(self.output_dir, 'circrna'),

        }

    def run(self):
        super(ExpMakecircTool, self).run()

        self.collect_count_matrix()
        self.set_output()
        self.end()



    def collect_count_matrix(self):
        if os.path.isdir(self.dir['count']):
            shutil.rmtree(self.dir['count'])
        os.makedirs(self.dir['count'])

        shutil.copy(self.option('c_count').path, os.path.join(self.dir['count'], 'C.reads.txt'))

    def set_output(self):
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)

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
        c_exp_df.to_csv(c_exp_table, sep='\t')

        c_stat_table = os.path.join(self.dir['count'], 'C.stat.txt')
        count = pd.read_table(self.option('c_count').path, index_col='seq_id')
        sample_id = dict()
        for sample in count.columns:
            sample_id[sample] = [x for x in count.index if count[sample][x] != 0]
            print len(sample_id[sample])
        detail = pd.read_table(self.option('c_detail').path, index_col='circrna_id')
        sample_host_gene_id = dict()
        for sample in count.columns:
            detail_none = detail.fillna("")
            host_gene_id = list(detail_none['host_gene_id'][sample_id[sample]])
            host_gene_id_new = ','.join(host_gene_id).split(',')
            host_gene_id_new_none = [x for x in host_gene_id_new if x != '']
            host_gene_id_new_none_set = list(set(host_gene_id_new_none))
            sample_host_gene_id[sample] = host_gene_id_new_none_set

        sample_circid_hostid = dict()
        circ_num = list()
        host_gene_num = list()
        totlecirc = list()
        totlehost = list()




        for sample in count.columns:
            circ_num.append(len(sample_id[sample]))
            host_gene_num.append(len(sample_host_gene_id[sample]))
            totlecirc.extend(sample_id[sample])
            totlecirc_duplicate = list(set(totlecirc))
            totlecirc_duplicate_num = len(totlecirc_duplicate)
            totlehost.extend(sample_host_gene_id[sample])
            totlehost_duplicate = list(set(totlehost))
            totlehost_duplicate_num = len(totlehost_duplicate)



        sample_circid_hostid['circrna'] = circ_num
        sample_circid_hostid['host_gene'] = host_gene_num
        sample_circid_hostid['sample'] = list(count.columns)
        print sample_circid_hostid
        stat = pd.DataFrame(sample_circid_hostid)
        totle = pd.DataFrame({'sample':'total','circrna':totlecirc_duplicate_num, 'host_gene':totlehost_duplicate_num}, index=[1])
        add_totle = stat.append(totle,ignore_index=True)
        add_totle.to_csv(c_stat_table, index=False, sep='\t')



        shutil.copytree(self.dir['count'], os.path.join(self.output_dir, 'count'))
        self.option('count_dir').set_path(os.path.join(self.output_dir, 'count'))




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
            'name': 'whole_transcriptome.formation.exp_makecirc',
            'instant': False,
            'options': {
                # 'exp_class': 'mlc',
                # 't_tpm': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/LargeGush/output/filter_by_express/classifyquant/T.tpm.txt',
                # 'g_tpm': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/LargeGush/output/filter_by_express/classifyquant/G.tpm.txt',
                # 't_fpkm': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/LargeGush/output/filter_by_express/classifyquant/T.fpkm.txt',
                # 'g_fpkm': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/LargeGush/output/filter_by_express/classifyquant/G.fpkm.txt',
                # 't_count': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/LargeGush/output/filter_by_express/classifyquant/T.count.txt',
                # 'g_count': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/LargeGush/output/filter_by_express/classifyquant/G.count.txt',
                # 'annot_table': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/Annotation/output/allannot_class/all_annot.xls',
                # 'relat_table': '/mnt/ilustre/users/sanger-dev/workspace/20191031/Longrna_tsg_36023/LargeGush/output/filter_by_express/filtered_file/trans_type.xls',
                'c_rpm': '/mnt/ilustre/users/sanger-dev/workspace/20191219/Circrna_tsg_36602/CircBrush/output/RPM.txt',
                'c_detail': '/mnt/ilustre/users/sanger-dev/workspace/20191219/Circrna_tsg_36602/CircBrush/output/detail.txt',
                'c_count': '/mnt/ilustre/users/sanger-dev/workspace/20191219/Circrna_tsg_36602/CircBrush/output/count.txt',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
