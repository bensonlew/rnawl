# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

import json
import os
import gevent
import gevent.pool
import datetime
import time  # 统计拷贝时间
from bson import ObjectId
from biocluster.api.database.base import Base
from gevent import Greenlet
from bson import CodecOptions, SON

class CopyDemo(Base):
    def __init__(self, old_task_id, new_task_id, new_project_sn, new_member_id, new_member_type):
        super(CopyDemo, self).__init__()
        self._project_type = 'metabolome'
        self._old_task_id = old_task_id
        self._new_task_id = new_task_id
        self._new_project_sn = new_project_sn
        self._new_member_id = new_member_id
        self._new_member_type = new_member_type
        self.group_id = ""
        self.diff_group_id = ""
        self.metab_table_id = ""
        self.diff_table_pos_id = ""
        self.diff_table_neg_id = ""
        self.diff_metab_set_id = ""
        self.org_metab_set_id = ""
        self.venn_metabset_id = []
        self.diff_table_id = ''

    def get_main_table(self):
        table_relation = self.db["table_relation"]
        main_tables = []
        results = table_relation.find({"is_detail":"n","is_workflow":"y","is_first":"n"})
        for each in results:
            if not each in main_tables:
                main_tables.append(each["table"])
        return main_tables

    def change_old_doc(self, collection, old_doc_dict):
        if 'is_demo' in old_doc_dict:
            old_doc_dict['is_demo'] = 1
        if 'task_id' in old_doc_dict:
            old_doc_dict['task_id'] = self._new_task_id
        if 'project_sn' in old_doc_dict:
            old_doc_dict['project_sn'] = self._new_project_sn
        if 'member_id' in old_doc_dict:
            old_doc_dict['member_id'] = self._new_member_id
        if 'group_id' in old_doc_dict:
            old_doc_dict['group_id'] = ObjectId(self.group_id)
        if 'metab_table_id' in old_doc_dict:
            old_doc_dict['metab_table_id'] = ObjectId(self.metab_table_id)
        if 'params' in old_doc_dict:
            params_str = old_doc_dict['params']
            if params_str and params_str != "null" and len(params_str) > 1:
                if isinstance(params_str, dict):
                    params_dict = params_str
                else:
                    params_dict = json.loads(params_str)
                if 'task_id' in params_dict.keys():
                    params_dict['task_id'] = self._new_task_id
                if 'diff_group_id' in params_dict.keys():
                    params_dict['diff_group_id'] = self.diff_group_id
                if 'metab_table' in params_dict.keys():   #20200928  拉取demo，params的metab_table 对应的是_id
                    params_dict['metab_table'] = self.metab_table_id
                if 'group_id' in params_dict.keys():
                    params_dict['group_id'] = self.group_id
                if collection == "metabset_venn":
                    self.venn_metabset_id = self.venn_metabset_id[0:6]
                    venn_metabset_ids = ",".join(self.venn_metabset_id)
                    params_dict['metabset'] = venn_metabset_ids
                else:
                    pass
                    #拉demo，参数反选中使用的 main_id 注释掉
                    # if 'metab_set' in params_dict.keys():
                    #     params_dict['metab_set'] = self.diff_metab_set_id
                    # if 'metabset' in params_dict.keys():
                    #     params_dict['metabset'] = self.diff_metab_set_id
                if 'diff_table' in params_dict.keys():
                    params_dict['diff_table'] = self.diff_table_id
                    # if params_dict['table_type'] == "pos":
                    #     params_dict['diff_table'] = self.diff_table_pos_id
                    # elif params_dict['table_type'] == "neg":
                    #     params_dict['diff_table'] = self.diff_table_neg_id
                params = json.dumps(params_dict, sort_keys=True, separators=(',', ':'))
                old_doc_dict['params'] = params
        return old_doc_dict

    def insert_copied_doc(self, table_name, new_doc_dict):
        collection = self.db[table_name]
        return collection.insert_one(new_doc_dict).inserted_id

    def get_target_doc_dict_list(self, collection):
        collection = self.db[collection]
        opts = CodecOptions(SON)
        collection_son = collection.with_options(codec_options=opts)
        result = collection_son.find({"task_id": self._old_task_id},{"_id":0})
        return result

    def copy_target_main_table(self, collection):
        doc_dict_list = self.get_target_doc_dict_list(collection)
        if not doc_dict_list:
            return
        for doc_dict in doc_dict_list:
            new_doc_dict = self.change_old_doc(collection, doc_dict)
            main_table_id = self.insert_copied_doc(collection, new_doc_dict)
            if collection == "specimen_group":
                self.group_id = str(main_table_id)
            if collection == "specimen_group_compare":
                self.diff_group_id = str(main_table_id)
            if collection == "exp":
                self.metab_table_id = str(main_table_id)
            if collection == "metab_set":
                if doc_dict["name"] == "DiffSet_mix":
                    self.diff_metab_set_id = str(main_table_id)
                elif doc_dict["name"] == "Set_Origin" or doc_dict["name"] == "Set_Raw":
                    self.org_metab_set_id =  str(main_table_id)
                else:
                    self.venn_metabset_id.append(str(main_table_id))
            if collection == "exp_diff":
                if 'params' in doc_dict:
                    params_str = doc_dict['params']
                    if isinstance(params_str, dict):
                        params_dict = params_str
                    else:
                        params_dict = json.loads(params_str)

                    self.diff_table_id = str(main_table_id)
                    # if params_dict['table_type'] == "pos":
                    #     self.diff_table_pos_id = str(main_table_id)
                    # elif params_dict['table_type'] == "neg":
                    #     self.diff_table_neg_id = str(main_table_id)

        print "{} ok".format(collection)

    def run(self):
        target_collections = self.get_main_table()
        if not target_collections:
            return
        print('we will copy the following main tables:', target_collections)
        self.copy_target_main_table('specimen_group')
        self.copy_target_main_table('specimen_group_compare')
        self.copy_target_main_table('exp')
        self.copy_target_main_table('exp_diff')
        self.copy_target_main_table('sg_task')
        self.copy_target_main_table('metab_set')
        pool = gevent.pool.Pool(100)
        pool.map(self.copy_target_main_table, target_collections)

if __name__ == '__main__':
    start_time = time.time()
    old_task_id = 'metabolome'
    new_task_id = 'copy_demo_test'
    new_project_sn = 'copy_demo_test'
    new_member_id = 'copy_demo_test'
    copy_demo = CopyDemoMongo(old_task_id, new_task_id, new_project_sn, new_member_id)
    copy_demo.run()
    end_time = time.time()
    print("total time: {}".format(end_time - start_time))
