# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mbio.api.database.small_rna.api_base import ApiBase
from biocluster.api.database.base import report_check
import os
import datetime
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
from bson.objectid import ObjectId
import numpy as np
import unittest

class MirnaEdit(ApiBase):
    def __init__(self, bind_object):
        super(MirnaEdit, self).__init__(bind_object)

    @report_check
    def add_mirna_edit(self, binomial_dir, task_id='small_rna', project_sn='small_rna', params=None, is_novel=False):

        # check the values of arguments
        if not os.path.isdir(binomial_dir):
            self.bind_object.set_error('{} is not a dir, abord'.format(binomial_dir))

        # check file and make sample_list
        info_file = os.path.join(binomial_dir, 'result.info.txt')
        sample_list = list()
        if not os.path.isfile(info_file):
            self.bind_object.set_error('can not find result.info.txt in {}, abord'.format(binomial_dir))
        else:
            with open(info_file) as f:
                for line in f:
                    binomial_xls, sample = line.strip().split('\t')
                    sample_list.append(sample)
                    if not os.path.isfile(binomial_xls):
                        self.bind_object.set_error('can not find {} in {}, abord'.format(binomial_xls, binomial_dir))

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

        # prepare mature_sites_dict
        sg_mirna_edit_sites = self.db['sg_mirna_edit_sites']
        mature_sites_list = list(sg_mirna_edit_sites.find({}, {'_id': 0}))
        hairpin_dict = dict()
        for each in mature_sites_list:
            if each['hairpin_name'] in hairpin_dict:
                hairpin_dict[each['hairpin_name']].update({each['mature_name']: each['mature_sites']})
            else:
                hairpin_dict[each['hairpin_name']] = {each['mature_name']: each['mature_sites']}

        # prepare detail_pd
        columns = ['miRNA_name', 'pre_miRNA_name', 'W', 'ME']
        for sample in sample_list:
            columns.extend([
                '{}_SP'.format(sample), # ME Position in matrue
                '{}_PP'.format(sample), # ME Position in hairpin
                '{}_MN'.format(sample), # Mismatch Number
                '{}_TN'.format(sample), # Total Number
                '{}_PV'.format(sample), # Row P-value
                '{}_PA'.format(sample) # Adjusted P-value
            ])
        detail_pd = pd.DataFrame([[None] * len(columns)], columns=columns)
        detail_pd = detail_pd.set_index('miRNA_name')
        bak_row = detail_pd.copy()

        # read info_file and complete detail_pd
        rf = open(info_file)
        for line in rf:
            binomial_xls, sample = line.strip().split('\t')
            matrix = pd.read_table(binomial_xls)
            matrix.set_index('miRNA_name', inplace=True)
            for h in matrix.index:
                if h not in hairpin_dict:
                    continue
                for m in hairpin_dict[h]:
                    if matrix.loc[h]['Location_inside_pre_miRNA'] in hairpin_dict[h][m]:
                        if m in detail_pd.index:
                            if isinstance(detail_pd.loc[m], pd.Series):
                                tmp_row = detail_pd.reindex([m])
                                if tmp_row.loc[m, 'pre_miRNA_name'] == h.lower() and tmp_row.loc[m, 'W'] == matrix.loc[h]['Mismatch_type'][0] and tmp_row.loc[m, 'ME'] == matrix.loc[h]['Mismatch_type'][-1]:
                                    tmp_row.loc[m, '{}_SP'.format(sample)] = hairpin_dict[h][m].index(matrix.loc[h]['Location_inside_pre_miRNA']) + 1
                                    tmp_row.loc[m, '{}_PP'.format(sample)] = matrix.loc[h]['Location_inside_pre_miRNA']
                                    tmp_row.loc[m, '{}_MN'.format(sample)] = matrix.loc[h]['Number_of_reads_with_the_mismatch']
                                    tmp_row.loc[m, '{}_TN'.format(sample)] = matrix.loc[h]['Total_number_of_reads_in_this_position']
                                    tmp_row.loc[m, '{}_PV'.format(sample)] = matrix.loc[h]['Raw_P_value']
                                    tmp_row.loc[m, '{}_PA'.format(sample)] = matrix.loc[h]['Bonferroni_P_value']
                                    detail_pd.loc[m] = tmp_row.loc[m]
                                else:
                                    tmp_row = bak_row.copy()
                                    tmp_row = tmp_row.reindex([m])
                                    tmp_row.loc[m, 'pre_miRNA_name'] = h.lower()
                                    tmp_row.loc[m, 'W'] = matrix.loc[h]['Mismatch_type'][0]
                                    tmp_row.loc[m, 'ME'] = matrix.loc[h]['Mismatch_type'][-1]
                                    tmp_row.loc[m, '{}_SP'.format(sample)] = hairpin_dict[h][m].index(matrix.loc[h]['Location_inside_pre_miRNA']) + 1
                                    tmp_row.loc[m, '{}_PP'.format(sample)] = matrix.loc[h]['Location_inside_pre_miRNA']
                                    tmp_row.loc[m, '{}_MN'.format(sample)] = matrix.loc[h]['Number_of_reads_with_the_mismatch']
                                    tmp_row.loc[m, '{}_TN'.format(sample)] = matrix.loc[h]['Total_number_of_reads_in_this_position']
                                    tmp_row.loc[m, '{}_PV'.format(sample)] = matrix.loc[h]['Raw_P_value']
                                    tmp_row.loc[m, '{}_PA'.format(sample)] = matrix.loc[h]['Bonferroni_P_value']
                                    detail_pd.append(tmp_row)
                            elif isinstance(detail_pd.loc[m], pd.DataFrame):
                                tmp_pd = detail_pd.loc[m].copy()
                                tmp_count = 0
                                for i in range(tmp_pd.shape[0]):
                                    if tmp_pd.iloc[i]['pre_miRNA_name'] == h.lower() and tmp_pd.iloc[i]['W'] == matrix.loc[h]['Mismatch_type'][0] and tmp_pd.iloc[i]['ME'] == matrix.loc[h]['Mismatch_type'][-1]:
                                        tmp_serie = tmp_pd.iloc[i]
                                        tmp_serie['{}_SP'.format(sample)] = hairpin_dict[h][m].index(matrix.loc[h]['Location_inside_pre_miRNA']) + 1
                                        tmp_serie['{}_SP'.format(sample)] = matrix.loc[h]['Location_inside_pre_miRNA']
                                        tmp_serie['{}_SP'.format(sample)] = matrix.loc[h]['Number_of_reads_with_the_mismatch']
                                        tmp_serie['{}_SP'.format(sample)] = matrix.loc[h]['Total_number_of_reads_in_this_position']
                                        tmp_serie['{}_SP'.format(sample)] = matrix.loc[h]['Raw_P_value']
                                        tmp_serie['{}_SP'.format(sample)] = matrix.loc[h]['Bonferroni_P_value']
                                        tmp_pd.iloc[i] = tmp_serie
                                        tmp_count += 1
                                if tmp_count > 1:
                                    self.bind_object.logger.critical('find situation problem in {}'.format(tmp_pd))
                                    debug_xls = os.path.join(self.bind_object.work_dir, 'debug.{}.xls'.format(datetime.datetime.now().strftime('%Y%m%d_%H%M%S')))
                                    detail_pd.to_csv(debug_xls, sep='\t')
                                    self.bind_object.logger.debug('generate debug file in {}'.format(debug_xls))
                                elif tmp_count == 1:
                                    detail_pd = detail_pd.drop(m)
                                    detail_pd.append(tmp_pd)
                                elif tmp_count == 0:
                                    tmp_row = bak_row.copy()
                                    tmp_row = tmp_row.reindex([m])
                                    tmp_row.loc[m, 'pre_miRNA_name'] = h.lower()
                                    tmp_row.loc[m, 'W'] = matrix.loc[h]['Mismatch_type'][0]
                                    tmp_row.loc[m, 'ME'] = matrix.loc[h]['Mismatch_type'][-1]
                                    tmp_row.loc[m, '{}_SP'.format(sample)] = hairpin_dict[h][m].index(matrix.loc[h]['Location_inside_pre_miRNA']) + 1
                                    tmp_row.loc[m, '{}_PP'.format(sample)] = matrix.loc[h]['Location_inside_pre_miRNA']
                                    tmp_row.loc[m, '{}_MN'.format(sample)] = matrix.loc[h]['Number_of_reads_with_the_mismatch']
                                    tmp_row.loc[m, '{}_TN'.format(sample)] = matrix.loc[h]['Total_number_of_reads_in_this_position']
                                    tmp_row.loc[m, '{}_PV'.format(sample)] = matrix.loc[h]['Raw_P_value']
                                    tmp_row.loc[m, '{}_PA'.format(sample)] = matrix.loc[h]['Bonferroni_P_value']
                                    detail_pd.append(tmp_row)
                        else:
                            tmp_row = bak_row.copy()
                            tmp_row = tmp_row.reindex([m])
                            tmp_row.loc[m, 'pre_miRNA_name'] = h.lower()
                            tmp_row.loc[m, 'W'] = matrix.loc[h]['Mismatch_type'][0]
                            tmp_row.loc[m, 'ME'] = matrix.loc[h]['Mismatch_type'][-1]
                            tmp_row.loc[m, '{}_SP'.format(sample)] = hairpin_dict[h][m].index(matrix.loc[h]['Location_inside_pre_miRNA']) + 1
                            tmp_row.loc[m, '{}_PP'.format(sample)] = matrix.loc[h]['Location_inside_pre_miRNA']
                            tmp_row.loc[m, '{}_MN'.format(sample)] = matrix.loc[h]['Number_of_reads_with_the_mismatch']
                            tmp_row.loc[m, '{}_TN'.format(sample)] = matrix.loc[h]['Total_number_of_reads_in_this_position']
                            tmp_row.loc[m, '{}_PV'.format(sample)] = matrix.loc[h]['Raw_P_value']
                            tmp_row.loc[m, '{}_PA'.format(sample)] = matrix.loc[h]['Bonferroni_P_value']
                            detail_pd = detail_pd.append(tmp_row)
        rf.close()
        detail_pd = detail_pd.drop([None])

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
            detail_pd.loc[m, 'edit_type'] = '{}to{}'.format(detail_pd.loc[m, 'W'], detail_pd.loc[m, 'ME'])
            wme = '{}{}'.format(detail_pd.loc[m, 'W'], detail_pd.loc[m, 'ME'])
            for sample in sample_list:
                if pd.isnull(detail_pd.loc[m, '{}_MN'.format(sample)]):
                    type_pd.loc[sample, wme] += 0
                else:
                    type_pd.loc[sample, wme] += detail_pd.loc[m, '{}_MN'.format(sample)]
        type_pd = type_pd.apply(np.int32)

        # create a new document in sg_mirna_edit_detail
        detail_pd = detail_pd.reset_index()
        detail_dict_list = detail_pd.to_dict('records')
        self.create_db_table('sg_mirna_edit_detail', detail_dict_list, tag_dict={'edit_id': main_id})

        # create a new document in sg_mirna_edit_type
        type_pd = type_pd.reset_index()
        type_dict_list = type_pd.to_dict('records')
        self.create_db_table('sg_mirna_edit_type', type_dict_list, tag_dict={'edit_id': main_id})

        # update sg_mirna_edit
        self.update_db_record('sg_mirna_edit', record_id=main_id, status='end')

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

        binomial_dir = '/mnt/ilustre/users/sanger-dev/workspace/20181128/Single_mirna_edit_7777_3070/MirnaEdit/output'
        wf.test_api.add_mirna_edit(binomial_dir)

if __name__ == '__main__':
    unittest.main()