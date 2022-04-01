# -*- coding: utf-8 -*-
# zouguanqing

import json
import os
import gevent
from gevent import pool
import datetime
import time  # 统计拷贝时间
from bson import ObjectId
from biocluster.api.database.base import Base
from gevent import Greenlet
from bson import CodecOptions, SON

class CopyDemo(Base):
    def __init__(self, old_task_id, new_task_id, new_project_sn, new_member_id, new_member_type):
        super(CopyDemo, self).__init__()
        self._project_type = 'meta'
        self._old_task_id = old_task_id
        self._new_task_id = new_task_id
        self._new_project_sn = new_project_sn
        self._new_member_id = new_member_id
        self._new_member_type = new_member_type
        self.sample_map = dict()
        self.group_map = dict()
        self.pair_group_map = dict()
        self.env_map = dict()
        self.env_group_map = dict()
        self.otu_id_map = dict()



    # 获取主表名称列表
    def get_main_table(self):
        table_relation = self.db["table_relation"]
        main_tables = []
        results = table_relation.find({"is_detail":"n"})
        for each in results:
            if not each in main_tables:
                main_tables.append(each["table"])
        return main_tables

    # 更新单个主表的数据
    def _change_old_doc(self, collection, old_doc_dict):
        if 'main_id' not in old_doc_dict:
            old_doc_dict['main_id'] = old_doc_dict['_id']
        ori_id = old_doc_dict.pop('_id')

        if 'is_demo' in old_doc_dict:
            old_doc_dict['is_demo'] = 1
        if 'task_id' in old_doc_dict:
            old_doc_dict['task_id'] = self._new_task_id
        if 'project_sn' in old_doc_dict:
            old_doc_dict['project_sn'] = self._new_project_sn
        if 'member_id' in old_doc_dict:
            old_doc_dict['member_id'] = self._new_member_id
        if 'group_id' in old_doc_dict:
            if old_doc_dict['group_id'] in ['All','all']:
                old_doc_dict['group_id'] = old_doc_dict['group_id']
            else:
                old_doc_dict['group_id'] = ObjectId(self.group_map[str(old_doc_dict['group_id'])])

        if 'otu_id' in old_doc_dict  and collection != 'sg_otu' :
            old_otu_id = str(old_doc_dict['otu_id'])
            old_doc_dict['otu_id'] = ObjectId(self.otu_id_map[old_otu_id])

        if 'params' in old_doc_dict:
            params_str = old_doc_dict['params']
            if params_str and params_str != "null" and len(params_str) > 1:
                if isinstance(params_str, dict):
                    params_dict = params_str
                else:
                    params_dict = json.loads(params_str)
                if 'task_id' in params_dict.keys():
                    params_dict['task_id'] = self._new_task_id

                if 'otu_id' in params_dict and collection != 'sg_otu':  # 后面再单独更新sg_otu 表中的params 的otu_id
                    old_otu_id = params_dict['otu_id']
                    params_dict['otu_id'] = self.otu_id_map[old_otu_id]

                if 'group_detail' in params_dict:
                    for group_name in params_dict['group_detail']:
                        for index in range(len(params_dict['group_detail'][group_name])):
                            old_name = params_dict['group_detail'][group_name][index]
                            new_name = self.sample_map[old_name]
                            params_dict['group_detail'][group_name][index] = new_name
                if 'env_id' in params_dict:
                    old_env_id = params_dict['env_id']
                    params_dict['env_id'] = self.env_map[old_env_id]

                if 'env_group_id' in params_dict:
                    old_env_group_id = params_dict['env_group_id']
                    params_dict['env_group_id'] = self.env_group_map[old_env_group_id]

                if 'group_id' in params_dict:
                    old_group_id = params_dict['group_id']
                    if old_group_id in ['All','all']:
                        params_dict['group_id'] = old_group_id
                    else:
                        params_dict['group_id'] = self.group_map[old_group_id]

                params = json.dumps(params_dict, sort_keys=True, separators=(',', ':'))
                old_doc_dict['params'] = params

        if collection == 'sg_specimen_group':
            new_specimen_names = []
            if not self.sample_map:
                raise Exception('需要先copy sg_specimen表')
            for group_info in old_doc_dict['specimen_names']:
                tmp = dict()
                for s in group_info:
                    if s not in self.sample_map:
                        raise Exception('%s 没有在sg_specimen中'%s)
                    new_s = self.sample_map[s]
                    tmp[new_s] = group_info[s]
                new_specimen_names.append(tmp)
            old_doc_dict['specimen_names'] = new_specimen_names

        elif collection == 'sg_specimen_pair_group':
            for type in ['specimen_list', 'pair_list']:
                new_sample_list = []
                for sample in old_doc_dict[type].split(','):
                    new_sample = self.sample_map[sample]
                    new_sample_list.append(new_sample)
                old_doc_dict[type] = ','.join(new_sample_list)
        elif collection == 'sg_task':
            if 'member_type' in old_doc_dict:
                old_doc_dict['member_type'] = self._new_member_type

        elif collection == 'sg_env_group':
            old_env_id = str(old_doc_dict['env_id'])
            new_env_id = self.env_map[old_env_id]
            old_doc_dict['env_id'] = ObjectId(new_env_id)

        elif collection == 'sg_status':
            if old_doc_dict['type_name'] == 'sg_otu':
                old_table_id = str(old_doc_dict['table_id'])
                old_doc_dict['table_id'] = ObjectId(self.otu_id_map[old_table_id])

        return old_doc_dict


    def _insert_copied_doc(self, table_name, new_doc_dict):
        collection = self.db[table_name]
        return collection.insert_one(new_doc_dict).inserted_id

    def get_target_doc_dict_list(self, collection, origin_filter=False):
        collection = self.db[collection]
        opts = CodecOptions(SON)
        collection_son = collection.with_options(codec_options=opts)

        if origin_filter:
            result = collection_son.find({"task_id": self._old_task_id, "name": {"$regex": "_Origin$"}})
        else:
            result = collection_son.find({"task_id": self._old_task_id})
        return result

    def copy_target_main_table(self, collection,update_status=False, origin_filter=True):   # 不要copy sg_status
        doc_dict_list = self.get_target_doc_dict_list(collection, origin_filter=origin_filter)
        if not doc_dict_list:
            return
        for doc_dict in doc_dict_list:
            new_doc_dict = self._change_old_doc(collection, doc_dict)
            main_table_id = self._insert_copied_doc(collection, new_doc_dict)
            print('insert %s : %s'%(collection, str(main_table_id)))

            if update_status:
                status_collection = self.db['sg_status']
                status_all = status_collection.find({"table_id": new_doc_dict['main_id'], "task_id": self._old_task_id})
                if status_all:
                    for status in status_all:
                        # table_id 等于 main_id 所以注释下2行
                        # if 'table_id' in status:
                        #     status['table_id'] = main_table_id
                        new_status_dict = self._change_old_doc('sg_status', status)
                        new_status_id = self._insert_copied_doc('sg_status', new_status_dict)
                        print("insert %s : %s"%('sg_status', str(new_status_id)))

            if collection == "sg_specimen":
                old_sample = str(new_doc_dict['main_id'])
                self.sample_map[old_sample] = str(main_table_id)
                self.db['sg_specimen'].update({"_id":main_table_id},{"$set":{"main_id":main_table_id}})
                self.copy_detail(old_sample, main_table_id, 'sg_specimen_step', 'specimen_id')
            elif collection == 'sg_specimen_group':
                old_group_id = str(new_doc_dict['main_id'])
                self.group_map[old_group_id] = str(main_table_id)
                self.db['sg_specimen_group'].update({"_id":main_table_id},{"$set":{"main_id":main_table_id}})

            elif collection =='sg_specimen_pair_group':
                old_pair_group_id = str(new_doc_dict['main_id'])
                self.pair_group_map[old_pair_group_id] = str(main_table_id)
                self.db['sg_specimen_pair_group'].update({"_id":main_table_id},{"$set":{"main_id":main_table_id}})

            elif collection == 'sg_env':
                old_env_id = str(new_doc_dict['main_id'])
                self.env_map[old_env_id] = str(main_table_id)
                details = self.db['sg_env_detail'].find({"env_id":new_doc_dict['main_id']},{"_id":0})
                insert_data = []
                for detail in details:
                    detail['env_id'] = main_table_id
                    insert_data.append(detail)
                self.db['sg_evn_detail'].insert_many(insert_data)
                self.db['sg_env'].update({"_id":main_table_id},{"$set":{"main_id":main_table_id}})

            elif collection == 'sg_env_group':
                old_env_group_id = str(new_doc_dict['main_id'])
                self.env_group_map[old_env_group_id] = str(main_table_id)
                self.db['sg_env_group'].update({"_id":main_table_id},{"$set":{"main_id":main_table_id}})


            elif collection == 'sg_otu':
                old_main_id = new_doc_dict['main_id']

                details = self.db['sg_otu_specimen'].find({"otu_id":old_main_id},{"_id":0})
                insert_data = []
                for detail in details:
                    detail['otu_id'] = main_table_id
                    old_sample_id = str(detail['specimen_id'])
                    detail['specimen_id'] = ObjectId(self.sample_map[old_sample_id])
                    insert_data.append(detail)
                self.db['sg_otu_specimen'].insert_many(insert_data)

                self.copy_detail(old_main_id, main_table_id, 'sg_otu_detail','otu_id')
                self.copy_detail(old_main_id, main_table_id, 'sg_otu_detail_level','otu_id')
                self.copy_detail(old_main_id, main_table_id, 'sg_otu_seq','otu_id')
                self.copy_detail(old_main_id, main_table_id, 'sg_otu_summary','otu_id')

                self.otu_id_map[str(old_main_id)] = str(main_table_id)
                self.db['sg_otu'].update({"_id":main_table_id},{"$set": {"main_id": main_table_id}})




        print "{} ok".format(collection)

    def copy_sg_task_workflow_params(self):
        ori_task_id_num = self._old_task_id.split('_')[-1]
        all_finds = self.db['sg_task_workflow_params'].find({"task_id": ori_task_id_num})
        new_task_id_num = self._new_project_sn.split('_')[-1]
        for table in all_finds:
            table['task_id'] = new_task_id_num
            table.pop('_id')
            self.db['sg_task_workflow_params'].insert_one(table)

    def copy_detail(self, ori_main_id, new_main_id,detail_collection, main_id_name):
        if not isinstance(ori_main_id, ObjectId):
            ori_main_id = ObjectId(ori_main_id)
        if not isinstance(new_main_id, ObjectId):
            new_main_id = ObjectId(new_main_id)

        detail_c = self.db[detail_collection]
        details = detail_c.find({main_id_name: ori_main_id})
        insert_datas = []
        for detail in details:
            detail.pop('_id')
            detail[main_id_name] = new_main_id
            if detail_collection in ['sg_specimen_step','sg_otu_detail_level']:
                detail['task_id'] = self._new_task_id
                detail['project_sn'] = self._new_project_sn
            insert_datas.append(detail)

        if insert_datas:
            detail_c.insert_many(insert_datas)


    def updata_sg_otu_otu_id(self):
        sg_otus = self.db['sg_otu'].find({"task_id" : self._new_task_id})
        for sg_otu in sg_otus:
            update_dic = {}
            if 'params' in sg_otu:
                params = json.loads(sg_otu['params'])
                old_otu_id = params['otu_id']
                params['otu_id'] = self.otu_id_map[old_otu_id]
                update_dic['params'] = json.dumps(params, sort_keys=True, separators=(',', ':'))
            if 'from_id' in sg_otu:
                new_from_id = self.otu_id_map[str(sg_otu['from_id'])]
                update_dic['from_id'] = ObjectId(new_from_id)
            if 'otu_id' in sg_otu:
                new_otu_id = self.otu_id_map[str(sg_otu['otu_id'])]
                update_dic['otu_id'] = ObjectId(new_otu_id)

            if update_dic:
                self.db['sg_otu'].update({"_id": sg_otu['_id']},{"$set": update_dic})


    def run(self):
        target_collections = self.get_main_table()
        if not target_collections:
            return

        # self.copy_sg_task_workflow_params()

        for c in ['sg_task_workflow_params']:
            if c in target_collections:
                target_collections.remove(c)

        for c in ['sg_task', 'sg_specimen','sg_specimen_group','sg_specimen_pair_group','sg_env',
                  'sg_env_group','sg_valid_sequence_info','sg_specimen_original_sequence',
                  'sg_raw_sequence_info','sg_otu']:  ##先复制样本，分组，在 sg_otu ,再是其他的主表
            if c == 'sg_otu':
                self.copy_target_main_table(c)
            else:
                self.copy_target_main_table(c, origin_filter=False)
            if c in target_collections:
                target_collections.remove(c)

        print('start updata_sg_otu_otu_id')
        self.updata_sg_otu_otu_id()

        print('we will copy the following main tables (*_Origin table) :', target_collections)

        p = pool.Pool(100)

        p.map(self.copy_target_main_table, target_collections)

if __name__ == '__main__':
    start_time = time.time()
    old_task_id = 'tsg_35938'
    new_task_id = 'copy_demo_tsg_35938'
    new_project_sn = 'copy_demo_test'
    new_member_id = 'copy_demo_test'

    copy_demo = CopyDemo(old_task_id, new_task_id, new_project_sn, new_member_id,'2')
    copy_demo.run()
    end_time = time.time()
    print("total time: {}".format(end_time - start_time))
