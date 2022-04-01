# -*- coding: utf-8 -*-

from bson.objectid import ObjectId
import datetime
from types import StringTypes
import os
import json
import pickle
import unittest
import pandas as pd
from biocluster.api.database.base import report_check
from mbio.api.database.whole_transcriptome.api_base import ApiBase
from bson.son import SON

class Synteny(ApiBase):
    def __init__(self, bind_object):
        super(Synteny, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_syteney_main(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
        if main_id is None:
            name = "synteny"+'_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v1",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='synteny',
                params=[],
                status="start",
            )
            main_id = self.create_db_table('synteny', [main_info])
        else:
            main_id = ObjectId(main_id)
        return main_id

    def add_detail(self, filepath, table, main_id, main_table,
                   columns=None, mongo_keys=None, header=True, header_lower=False,
                   tag_key=[], tag_value=[], type=None):
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)
        h = 0 if header else None
        df = pd.read_csv(filepath, header=h, index_col=False, sep='\t')
        if columns:
            df = df[columns]
        if mongo_keys:
            df = df.rename(columns=mongo_keys)
        df[main_table] = main_id
        if type:
            df["type"] = type
        if tag_key:
            if len(tag_key) != len(tag_value):
                self.set_error('新加的tag key列表必须和value列表等长')
            for i in range(len(tag_key)):
                df[tag_key[i]] = tag_value[i]
        if header_lower:
            df.columns = map(lambda x: x.lower(), df.columns)
        data = self.df_to_mongo(df)
        self._run_add(table, data)

    def df_to_mongo(self, df):
        keys = map(lambda x: x.replace('.', '_'), df.columns)
        mongo_data = []
        df.apply(
            lambda x: mongo_data.append(SON(dict(zip(keys, x)))), axis=1
        )
        return mongo_data

    def update_sg_status(self, table_id, data):
        self.update_table('sg_status', table_id, data, search_id='table_id')

    def update_table(self, table, table_id, data, search_id='_id'):
        if not isinstance(table_id, ObjectId):
            if isinstance(table_id, str):
                table_id = ObjectId(table_id)
            else:
                self.bind_object.set_error("{}必须为ObjectID对象或其字符串形式".format(table_id))

        tb = self.db[table]
        tb.update({search_id: table_id},
                  {'$set': data})

    def _run_add(self, name, data, main=False):
        try:
            if main:
                collection = self.db[name]
                main_id = collection.insert_one(SON(data)).insert_id
                collection.update_one({'_id': main_id},
                                      {'$set': {'main_id': main_id}},
                                      )
                return main_id
            else:
                collection = self.db[name]
                collection.insert_many(data)
        except Exception, e:
            self.bind_object.set_error('导入表%s出错：%s') % (name, e)
        else:
            self.bind_object.logger.info('导入表%s成功！' % name)

    def updata_status(self,main_id):
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)
        dotplot_data_dict = {"query_super_id":"query_super_id","query":"query","query_start":"query_start","query_end":"query_end",
                             "ref_super_id":"ref_super_id","ref":"ref","ref_start":"ref_start","ref_end":"ref_end",
                             "identity":"identity","direction":"direction","query_scaf":"query_scaf","ref_scaf":"ref_scaf","condition":{"type":"synteny"}}
        dotplot_data_info = json.dumps(dotplot_data_dict, sort_keys=False, separators=(',', ':'))
        synteny_data_dict = {"query_len":"query_len","query_end":"query_end","query_scaf":"query_scaf",
                             "query_super_end":"query_super_end","ref_super_end":"ref_super_end","query_start":"query_start",
                             "ref_start":"ref_start","ref_super_start":"ref_super_start","query_super_id":"query_super_id",
                             "ref_scaf":"ref_scaf","query_super_start":"query_super_start","ref_end":"ref_end","ref_len":
                                 "ref_len","query":"query","ref":"ref","identity":"identity","ref_super_id":"ref_super_id","condition":{"type":"synteny"}}
        synteny_data_info = json.dumps(synteny_data_dict, sort_keys=False, separators=(',', ':'))
        band_data_dict = {"sp":"sp","chr":"chr","len":"len","condition":{"type":"band_synteny"}}
        band_data_info = json.dumps(band_data_dict, sort_keys=False, separators=(',', ':'))
        try:
            self.update_db_record('synteny', main_id, status="end", main_id=main_id,
                                 synteny_data=synteny_data_info, band_data=band_data_info,dotplot_data=dotplot_data_info,)
        except Exception as e:
            self.bind_object.logger.error("导入synteny数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入synteny数据成功")








