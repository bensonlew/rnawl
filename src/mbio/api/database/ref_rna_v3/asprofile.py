# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
from bson.objectid import ObjectId
import datetime
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase


class Asprofile(ApiBase):
    def __init__(self, bind_object):
        super(Asprofile, self).__init__(bind_object)
        self._project_type = 'ref_rna_v2'

    def add_asprofile_result(self, result_table, s3_path, main_id=None, task_id='ASprofile', project_sn='ASprofile'):
        if main_id is None:
            name = "ASprofile"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            # if type(params) == dict:
            #     params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v3",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='ASprofile main table',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('sg_asprofile', [main_info])
        else:
            if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
                main_id = ObjectId(main_id)
        df = pd.read_table(result_table, header=0, sep='\t', keep_default_na=False)
        event_type = list(set(list(df['event_type'])))
        df['asprofile_id'] = main_id
        self.create_db_table('sg_asprofile_result_detail', df.to_dict('r'))
        self.update_db_record('sg_asprofile', main_id, status="end", event_type=event_type, asprofile_path=s3_path)
        return main_id
    def add_asprofile_statistics(self, statistics_table, main_id, sample_list, task_id='ASprofile', project_sn='ASprofile'):
        if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
            main_id = ObjectId(main_id)
        df = pd.read_table(statistics_table, header=0, sep='\t', keep_default_na=False)
        df_1 = df.set_index('sample')
        df_2 = df_1.reindex(sample_list.split(','))
        df_3 = df_2.reset_index()
        df_3['asprofile_id'] = main_id
        df_3_list = df_3.to_dict('r')
        df_3_list.sort(key=lambda x: sample_list.split(',').index(x['sample']))

        self.create_db_table('sg_asprofile_statistics_detail', df_3_list)


        self.update_db_record('sg_asprofile', main_id, status="end", main_id=main_id)

    def add_diff_asprofile_result(self, result_table, main_id=None, task_id='ASprofile', project_sn='ASprofile'):
        if main_id is None:
            name = "Diff_ASprofile"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            # if type(params) == dict:
            #     params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v3",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='Diff_ASprofile main table',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('sg_asprofile_diff', [main_info])
        else:
            if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
                main_id = ObjectId(main_id)
        df = pd.read_table(result_table, header=0, sep='\t', keep_default_na=False)
        event_type = list(set(list(df['event_type'])))
        df['diff_id'] = main_id
        self.create_db_table('sg_asprofile_diff_detail', df.to_dict('r'))
        self.update_db_record('sg_asprofile_diff', main_id, status="end", event_type=event_type)
        return main_id