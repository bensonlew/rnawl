# -*- coding: utf-8 -*-
from biocluster.api.database.base import Base, report_check
import pandas as pd
import datetime
import types
from bson.objectid import ObjectId
import json


class Metabset(Base):
    def __init__(self, bind_object):
        super(Metabset, self).__init__(bind_object)
        self._project_type = "metabolome"

    @report_check
    def add_metab_set(self, name, me_num,not_delete=False):
        # me_num 代谢物个数
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': "代谢集",
            'created_ts': created_ts,
            # "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            'status': 'end',
            'is_use': 2,  # 工作流穿件的代谢集不可删
            'metabset_length': me_num,
            'name': name
        }
        if not_delete:
            insert_data['not_delete'] = 1
        collection = self.db['metab_set']
        main_id = collection.insert_one(insert_data).inserted_id
        collection.update_one({'_id': main_id}, {'$set': {'main_id': main_id}})
        return main_id

    @report_check
    def add_metab_set_detail(self, main_id, metab_list):
        # metab_list 代谢物的集合列表
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="54700401")
        insert_data = {
            'set_id': main_id,
            'set_list': metab_list
        }
        collection = self.db['metab_set_detail']
        collection.insert_one(insert_data)

    def run(self):
        main_id = self.add_metab_set("Set_Origin", me_num)
        self.add_metab_set_detail(main_id, metab_list)
