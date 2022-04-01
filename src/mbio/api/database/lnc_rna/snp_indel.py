# -*- coding: utf-8 -*-
# __author__ = 'qindanhua, qinjincheng'

from mbio.api.database.ref_rna_v2.api_base import ApiBase
from biocluster.api.database.base import report_check
import datetime
import json
from bson.objectid import ObjectId
import os
import pandas as pd
import unittest

class SnpIndel(ApiBase):
    '''
    last_modify: 2019.04.23
    '''
    def __init__(self, bind_object):
        super(SnpIndel, self).__init__(bind_object)
        self._project_type = 'lnc_rna'

    @report_check
    def add_snp(self, module_output, method_type, upload_dir, gene_type_tsv, params=None, main_id=None):
        if main_id is None:
            project_sn = self.bind_object.sheet.project_sn
            task_id = self.bind_object.sheet.id
            time_now = datetime.datetime.now()
            created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
            name = 'Snp_{}_{}'.format(method_type, time_now.strftime('%Y%m%d_%H%M%S'))
            if type(params) == dict:
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = {
                'project_sn': project_sn,
                'task_id': task_id,
                'created_ts': created_ts,
                'name': name,
                'desc': 'Snp main table',
                'params': params,
                'status': 'start'
            }
            snp_id = self.create_db_table('sg_snp', [main_info])
        else:
            snp_id = ObjectId(main_id)
        ret = self.add_snp_detail(module_output, snp_id, upload_dir, gene_type_tsv)
        if ret:
            insert_dict = {'main_id': snp_id, 'status': 'end'}
            for line in open(os.path.join(module_output, 'main_info.txt')):
                items = line.strip().split('\t')
                if len(items) >= 2:
                    insert_dict.update({items[0]: items[1].split(';')})
            self.update_db_record('sg_snp', snp_id, insert_dict=insert_dict)
        return snp_id

    def add_snp_detail(self, module_output, snp_id, upload_dir, gene_type_tsv):
        gene_type_dict = dict()
        for line in open(gene_type_tsv):
            items = line.strip().split('\t')
            if len(items) > 3:
                gene_type_dict[items[0]] = items[2]
        data_anno = os.path.join(module_output, 'data_anno_pre.xls')
        with open(data_anno) as f:
            headers = f.readline().strip().split('\t')
            data_list = list()
            for row, line in enumerate(f):
                data = dict()
                line = line.strip().split('\t')
                if len(headers) == len(line):
                    for n in range(len(headers)):
                        if headers[n].endswith('_mut_rate') or headers[n] in ['qual']:
                            try:
                                data.update({headers[n]: float(line[n])})
                            except:
                                data.update({headers[n]: line[n]})
                        elif headers[n].endswith('_reads_rate') or headers[n] in ['start', 'end', 'reads_num']:
                            try:
                                data.update({headers[n]: int(line[n])})
                            except:
                                data.update({headers[n]: line[n]})
                        else:
                            data.update({headers[n]: line[n]})
                    data.update({'snp_id': snp_id})
                    if data['gene'] in gene_type_dict:
                        data.update({'rna_type': gene_type_dict[data['gene']]})
                    else:
                        data.update({'rna_type': 'other'})
                    data_list.append(data)
                    if row != 0 and row % 100000 == 0:
                        self.create_db_table('sg_snp_detail', data_list)
                        data_list = list()
            if data_list:
                self.create_db_table('sg_snp_detail', data_list)
                data_list = list()
        os.rename(os.path.join(module_output, 'data_anno_pre.xls'), os.path.join(module_output, 'snp_annotation.xls'))

        df_sds = pd.read_table(os.path.join(module_output, 'snp_depth_statistics.xls'))
        df_sds = df_sds.fillna('')
        df_sds['snp_id'] = snp_id
        depth_data = df_sds.to_dict('records')
        for depth in depth_data:
            for key in depth:
                try:
                    depth[key] = float(depth[key])
                except:
                    pass
        self.create_db_table('sg_snp_stat', depth_data)
        df_sds.rename(columns={'range_key': 'Depth'}, inplace=True)
        df_sds.drop(['snp_id', 'stat_type'], axis=1, inplace=True)
        df_sds_select_col = df_sds.columns[0:-1].tolist()
        df_sds_select_col.insert(0, df_sds.columns[-1])
        df_depth = df_sds[df_sds_select_col]
        snp_depth_statistics = os.path.join(upload_dir, 'snp_depth_statistics.xls')
        df_depth.to_csv(snp_depth_statistics, sep='\t', header=True, index=False)

        df_sfs = pd.read_table(os.path.join(module_output, 'snp_freq_statistics.xls'))
        df_sfs = df_sfs.fillna('')
        df_sfs['snp_id'] = snp_id
        freq_data = df_sfs.to_dict('records')
        for freq in freq_data:
            for key in freq:
                try:
                    freq[key] = float(freq[key])
                except:
                    pass
        self.create_db_table('sg_snp_stat', freq_data)
        df_sfs.rename(columns={'range_key': 'Freq'}, inplace=True)
        df_sfs.drop(['snp_id', 'stat_type'], axis=1, inplace=True)
        df_sfs_select_col = df_sfs.columns[0:-1].tolist()
        df_sfs_select_col.insert(0, df_sfs.columns[-1])
        df_freq_data = df_sfs[df_sfs_select_col]
        snp_freq_statistics = os.path.join(upload_dir, 'snp_freq_statistics.xls')
        df_freq_data.to_csv(snp_freq_statistics, sep='\t', header=True, index=False)

        df_stts = pd.read_table(os.path.join(module_output, 'snp_transition_tranversion_statistics.xls'))
        df_stts = df_stts.fillna('')
        df_stts['snp_id'] = snp_id
        type_data = df_stts.to_dict('records')
        for type in type_data:
            for key in type:
                try:
                    type[key] = float(type[key])
                except:
                    pass
        self.create_db_table('sg_snp_stat', type_data)
        df_stts.rename(columns={'range_key': 'Type'}, inplace=True)
        df_stts.drop(['snp_id', 'stat_type'], axis=1, inplace=True)
        df_stts_select_col = df_stts.columns[0:-1].tolist()
        df_stts_select_col.insert(0, df_stts.columns[-1])
        df_type_data = df_stts[df_stts_select_col]
        snp_transition_tranversion_statistics = os.path.join(upload_dir, 'snp_transition_tranversion_statistics.xls')
        df_type_data.to_csv(snp_transition_tranversion_statistics, sep='\t', header=True, index=False)

        df_spd = pd.read_table(os.path.join(module_output, 'snp_position_distribution.xls'))
        df_spd = df_spd.fillna('')
        df_spd['snp_id'] = snp_id
        snp_region_data = df_spd.to_dict('records')
        for pos in snp_region_data:
            for key in pos:
                try:
                    pos[key] = float(pos[key])
                except:
                    pass
        self.create_db_table('sg_snp_stat', snp_region_data)
        df_spd.rename(columns={'range_key': 'Region'}, inplace=True)
        df_spd.drop(['snp_id', 'stat_type'], axis=1, inplace=True)
        df_spd_select_col = df_spd.columns[0:-1].tolist()
        df_spd_select_col.insert(0, df_spd.columns[-1])
        df_snp_pos_data = df_spd[df_spd_select_col]
        snp_position_distribution = os.path.join(upload_dir, 'snp_position_distribution.xls')
        df_snp_pos_data.to_csv(snp_position_distribution, sep='\t', header=True, index=False)

        df_ipd = pd.read_table(os.path.join(module_output, 'indel_position_distribution.xls'))
        df_ipd = df_ipd.fillna('')
        df_ipd['snp_id'] = snp_id
        indel_region_data = df_ipd.to_dict('records')
        for indel in indel_region_data:
            for key in indel:
                try:
                    indel[key] = float(indel[key])
                except:
                    pass
        self.create_db_table('sg_snp_stat', indel_region_data)
        df_ipd.rename(columns={'range_key': 'Region'}, inplace=True)
        df_ipd.drop(['snp_id', 'stat_type'], axis=1, inplace=True)
        df_ipd_select_col = df_ipd.columns[0:-1].tolist()
        df_ipd_select_col.insert(0, df_ipd.columns[-1])
        df_indel_pos_data = df_ipd[df_ipd_select_col]
        indel_position_distribution = os.path.join(upload_dir, 'indel_position_distribution.xls')
        df_indel_pos_data.to_csv(indel_position_distribution, sep='\t', header=True, index=False)

        return True

class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.lnc_rna.lnc_rna_test_api import LncRnaTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'snp_indel_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'lnc_rna.lnc_rna_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = LncRnaTestApiWorkflow(wheet)
        wf.sheet.id = 'lnc_rna'
        wf.sheet.project_sn = 'lnc_rna'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('lnc_rna.snp_indel')
        params = {
            'method_type': 'samtools',
            'submit_location': 'snp',
            'task_id': 'lnc_rna',
            'task_type': 2,
        }
        wf.test_api.add_snp(
            module_output='/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/snp_indel/output',
            method_type='samtools',
            params=params,
            upload_dir='/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/snp_indel/upload',
        )

if __name__ == '__main__':
    unittest.main()