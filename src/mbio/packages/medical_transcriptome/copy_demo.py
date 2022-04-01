# -*- coding: utf-8 -*-
# __author__ = 'gdq'
import time
from biocluster.api.database.base import Base
import gevent.pool
import gevent.monkey
import json
from bson import CodecOptions, SON

#gevent.monkey.patch_all()


class CopyDemoMongo(Base):
    def __init__(self, old_task_id, new_task_id, new_project_sn, new_member_id, project_type='medical_transcriptome'):
        super(CopyDemoMongo, self).__init__()
        self._project_type = project_type
        self._old_task_id = old_task_id
        self._new_task_id = new_task_id
        self._new_project_sn = new_project_sn
        self._new_member_id = new_member_id

    def change_old_doc(self, old_doc_dict):
        if 'task_id' in old_doc_dict:
            old_doc_dict['task_id'] = self._new_task_id
        if 'project_sn' in old_doc_dict:
            old_doc_dict['project_sn'] = self._new_project_sn
        if 'member_id' in old_doc_dict:
            old_doc_dict['member_id'] = self._new_member_id
        if 'params' in old_doc_dict:
            params_str = old_doc_dict['params']
            if params_str and params_str != "null" and params_str != "none" and len(params_str) > 1:
                if isinstance(params_str, dict):
                    if 'task_id' in params_str:
                        old_doc_dict['params']['task_id'] = self._new_task_id
                else:
                    try:
                        params_dict = json.loads(params_str)
                        if 'task_id' in params_dict:
                            params_dict['task_id'] = self._new_task_id
                            params = json.dumps(params_dict, sort_keys=True, separators=(',', ':'))
                            old_doc_dict['params'] = params
                    except:
                        pass
        return old_doc_dict

    def insert_copied_doc(self, table_name, new_doc_dict):
        collection = self.db[table_name]
        return collection.insert_one(new_doc_dict).inserted_id

    def get_target_collection_names(self):
        """
        表名如果以'_detail'结尾，则认为是详情表，不予复制；表的记录字段中如果含有'task_id'则认为是主表，需复制！
        :return:
        """
        target_collections = []
        collections = self.db.collection_names()
        for coll in collections:
            # if coll.endswith('_detail') and coll != 'sg_t2g_detail':
            if coll.endswith('_detail'):
                continue
            if coll == "system.profile":
                continue
            doc_dict = self.db[coll].find_one()
            if doc_dict:
                if 'task_id' in doc_dict:
                    target_collections.append(coll)
        return target_collections

    def get_target_doc_dict_list(self, collection):
        collection = self.db[collection]
        opts = CodecOptions(SON)
        collection_son = collection.with_options(codec_options=opts)
        result = collection_son.find({"task_id": self._old_task_id}, {"_id": 0})
        return result

    def copy_target_main_table(self, collection):
        doc_dict_list = self.get_target_doc_dict_list(collection)
        if not doc_dict_list:
            return
        for doc_dict in doc_dict_list:
            new_doc_dict = self.change_old_doc(doc_dict)
            self.insert_copied_doc(collection, new_doc_dict)

    def run(self):
        target_collections = self.get_target_collection_names()
        if not target_collections:
            return
        print('we will copy the following main tables:', target_collections)
        pool = gevent.pool.Pool(100)
        pool.map(self.copy_target_main_table, target_collections)


if __name__ == '__main__':
    start_time = time.time()
    old_task_id = 'ref_rna_v2'
    new_task_id = 'copy_demo_test'
    new_project_sn = 'copy_demo_test'
    new_member_id = 'copy_demo_test'
    copy_demo = CopyDemoMongo(old_task_id, new_task_id, new_project_sn, new_member_id)
    copy_demo.run()
    end_time = time.time()
    print("total time: {}".format(end_time - start_time))

