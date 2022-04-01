# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
from bson.objectid import ObjectId
import datetime
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.medical_transcriptome.api_base import ApiBase


class Asprofile(ApiBase):
    def __init__(self, bind_object):
        super(Asprofile, self).__init__(bind_object)
        self._project_type = 'medical_transcriptome'

    def add_asprofile_result(self, result_table,s3_path, main_id=None, task_id='medical_transcriptome', project_sn='ASprofile',params = None ):
        if main_id is None:
            name = "ASprofile"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            # if type(params) == dict:
            #     params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='ASprofile main table',
                params = params,
                status = "start",
            )
            main_id = self.create_db_table('sg_asprofile', [main_info])
        else:
            if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
                main_id = ObjectId(main_id)
        # df = pd.read_table(result_table, header=0, sep='\t', keep_default_na=False)
        # event_type = list(set(list(df['event_type'])))

        splice_pd = pd.read_table(result_table, header=0, sep='\t')
        categories = list(splice_pd.columns[:-1])
        # df['asprofile_id'] = main_id
        # self.create_db_table('sg_asprofile_result_detail', df.to_dict('r'))
        self.update_db_record('sg_asprofile', main_id, status="end", asprofile_path=s3_path,event_type=categories)
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


    def add_asprofile_search_result(self,as_search_results, main_id=None, task_id='medical_transcriptome', project_sn='ASprofile_search'):
        if main_id is None:
            name = "ASprofile_search"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            # if type(params) == dict:
            #     params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='ASprofile search main table',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('sg_asprofile_search', [main_info])
        else:
            if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
                main_id = ObjectId(main_id)
        df = pd.read_table(as_search_results, header=0, sep='\t', keep_default_na=False)
        # event_type = list(set(list(df['event_type'])))
        df['asprofile_search_id'] = main_id
        #仅保留前50000行数据
        if df.shape[0] >= 50000:
            flag = 1
        else:
            flag =0
        if flag == 0:
            try:
                self.create_db_table('sg_asprofile_search_detail', df.to_dict('r'))
            except Exception as e:
                self.bind_object.set_error("导入cog表格：%s出错:%s" % (as_search_results, e))
            self.update_db_record('sg_asprofile_search', main_id, status="end")
        else:
            self.update_db_record('sg_asprofile_search', main_id, status="end",detail_show = "no")
        collection = self.db["sg_asprofile_search"]
        task_id = collection.find_one({"main_id": main_id, "status": "end"})["task_id"]
        num = collection.find({"task_id": task_id, "status": "end"}).count()
        #检查是否超过3条搜索记录(需先导入成功后删除)
        if num == 4:
            self.bind_object.logger.info("已有四条搜索记录,删除最早的搜索记录!" )
            deleted_record = collection.find_one({"task_id": task_id, "status": "end"})
            deleted_main_id = deleted_record["main_id"]
            try:
                collection.remove({"main_id": deleted_main_id, "task_id": task_id})
                detail_collection = self.db["sg_asprofile_search_detail"]
                detail_collection.delete_many({"asprofile_search_id": deleted_main_id})
            except:
                self.bind_object.set_error("删除记录失败")

        return main_id