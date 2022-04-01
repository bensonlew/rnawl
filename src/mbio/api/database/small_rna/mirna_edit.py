# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mbio.api.database.small_rna.api_base import ApiBase
from biocluster.api.database.base import report_check
import os
import datetime
import pandas as pd
import numpy as np
import unittest

class MirnaEdit(ApiBase):
    def __init__(self, bind_object):
        super(MirnaEdit, self).__init__(bind_object)

    @report_check
    def add_mirna_edit(self, result_dir, task_id='small_rna', project_sn='small_rna', params=None, is_novel=False):

        # check the values of arguments
        if not os.path.isdir(result_dir):
            self.bind_object.set_error('{} is not a dir, abord'.format(result_dir))

        # check file and make sample_list
        info_file = os.path.join(result_dir, 'result.info.txt')
        sample_list = list()
        if not os.path.isfile(info_file):
            self.bind_object.set_error('can not find result.info.txt in {}, abord'.format(result_dir))
        else:
            with open(info_file) as f:
                for line in f:
                    result_xls, sample = line.strip().split('\t')
                    sample_list.append(sample)
                    if not os.path.isfile(result_xls):
                        self.bind_object.set_error('can not find {} in {}, abord'.format(result_xls, result_dir))

        # prepare main table info
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        name = 'MirnaEdit_{}'.format(time_now.strftime('%Y%m%d_%H%M%S'))
        desc = 'mirna_edit_main_table'

        # insert one document to sg_mirna_edit, return _id as main_id
        main_info = {
            # essential keys
            'task_id': task_id,
            'project_sn': project_sn,
            'created_ts': created_ts,
            'name': name,
            'desc': desc,
            'params': params,
            # alternative keys
            'sample_list': sample_list,
            'is_novel': is_novel,
            # status of process
            'status': 'start'
        }
        main_id = self.create_db_table('sg_mirna_edit', [main_info])
        self.update_db_record('sg_mirna_edit', main_id, main_id=main_id)
        self.bind_object.logger.info('succeed in creating document in sg_mirna_edit')

        # read info_file and make detail_pd
        pd_dict = dict()
        for line in open(info_file):
            result_xls, sample = line.strip().split('\t')
            pd_dict[sample] = pd.read_table(result_xls)
            name_map_dict = {
                'mature_id': 'miRNA_name',
                'hairpin_id': 'pre_miRNA_name',
                'w': 'W',
                'e': 'E',
                'mp': 'MP',
                'hp': 'PP',
                'mis_num': '{}_MN'.format(sample),
                'total_num': '{}_TN'.format(sample),
                'p_value': '{}_PV'.format(sample),
                'p_adjust': '{}_PA'.format(sample)
            }
            pd_dict[sample].rename(columns=name_map_dict, inplace=True)
        if len(sample_list) >= 2:
            detail_pd = pd.merge(pd_dict[sample_list[0]], pd_dict[sample_list[1]], how='outer')
            for sample in sample_list[2:]:
                detail_pd = detail_pd.merge(pd_dict[sample], how='outer')
        else:
            detail_pd = pd_dict[sample_list[0]]
        all_result_pd = detail_pd.copy()

        # convert NaN to None
        for x in range(detail_pd.shape[0]):
            for y in range(detail_pd.shape[1]):
                if pd.isnull(detail_pd.iloc[x, y]):
                    detail_pd.iloc[x, y] = None

        # prepare type_pd
        type_row = sample_list
        type_col = ['AC', 'AG', 'AU', 'CA', 'CG', 'CU', 'GA', 'GC', 'GU', 'UA', 'UC', 'UG']
        type_arr = np.zeros((len(type_row), len(type_col)))
        type_pd = pd.DataFrame(type_arr, index=type_row, columns=type_col)
        type_pd.index.name = 'sample'

        # read detail_pd then complete detail_pd and type_pd
        for m in detail_pd.index:
            detail_pd.loc[m, 'edit_type'] = '{}to{}'.format(detail_pd.loc[m, 'W'], detail_pd.loc[m, 'E'])
            we = '{}{}'.format(detail_pd.loc[m, 'W'], detail_pd.loc[m, 'E'])
            for sample in sample_list:
                if pd.isnull(detail_pd.loc[m, '{}_MN'.format(sample)]):
                    type_pd.loc[sample, we] += 0
                else:
                    # type_pd.loc[sample, we] += detail_pd.loc[m, '{}_MN'.format(sample)]
                    type_pd.loc[sample, we] += 1
        type_pd = type_pd.apply(np.int32)

        # create a new document in sg_mirna_edit_detail
        detail_pd = detail_pd.reset_index()
        detail_dict_list = detail_pd.to_dict('records')

        # convert nan to None
        for i in range(len(detail_dict_list)):
            for k in detail_dict_list[i]:
                if pd.isnull(detail_dict_list[i][k]):
                    detail_dict_list[i][k] = '-'
        self.create_db_table('sg_mirna_edit_detail', detail_dict_list, tag_dict={'edit_id': main_id})

        # create a new document in sg_mirna_edit_type
        type_pd = type_pd.reset_index()
        type_dict_list = type_pd.to_dict('records')
        self.create_db_table('sg_mirna_edit_type', type_dict_list, tag_dict={'edit_id': main_id})

        # update sg_mirna_edit
        self.update_db_record('sg_mirna_edit', record_id=main_id, status='end')

        # create all.result.xls and remove result.info.txt
        os.remove(info_file)
        all_result_pd.to_csv(os.path.join(result_dir, 'all.result.xls'), sep='\t', index=None)

class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.small_rna.small_rna_test_api import SmallRnaTestApiWorkflow
        from biocluster.wsheet import Sheet

        data = {
            'id': 'mirna_edit_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'small_rna.small_rna_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = SmallRnaTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('small_rna.mirna_edit')

        result_dir = '/mnt/ilustre/users/sanger-dev/workspace/20200825/Smallrna_tsg_218510/MirnaEdit/output'
        wf.test_api.add_mirna_edit(result_dir)

if __name__ == '__main__':
    unittest.main()
