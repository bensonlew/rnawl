# -*- coding: utf-8 -*-
# __author__ = 'gdq'
import time
from biocluster.config import Config
import json
from bson import CodecOptions, SON
from concurrent.futures import ThreadPoolExecutor


class DeleteRecordsMongo(object):
    def __init__(self, task_id, project_type, submit_location=all, status='failed'):
        self.db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        self.task_id = task_id
        self.status = status
        self._project_type = project_type
        self.submit_location = submit_location
        self.delete_set = self.get_delete_targets()
        if submit_location == all:
            self.delete_targets = self.delete_set
        else:
            self.delete_targets = self.get_some_delete_target(submit_location)

    def remove_sg_status(self):
        conn = self.db['sg_status']
        a = conn.find({"task_id": self.task_id, "status": self.status, "submit_location": self.submit_location}, {'_id': 1})
        if type(a) == dict:
            a = [a]
        for each in a:
            self.remove_db_record('sg_status', _id=each['_id'])

    def get_some_delete_target(self, submit_location):
        delete_targets = []
        target_collection = ''
        collections = self.db.collection_names()
        for coll in collections:
            # if coll.endswith('_detail') and coll != 'sg_t2g_detail':
            if coll.endswith('_detail'):
                continue
            doc_dicts = self.db[coll].find().sort([('_id',-1)])
            doc_dict = list()
            try:
                doc_dict = list(doc_dicts)[0]
            except:
                pass
            if doc_dict:
                if 'params' in doc_dict:
                    params_str = doc_dict['params']
                    if params_str and params_str != "null" and len(params_str) > 1:
                        if isinstance(params_str, dict):
                            if 'submit_location' in params_str:
                                if doc_dict['params']['submit_location'] == submit_location:
                                    target_collection = coll
                                    break
                        else:
                            try:
                                params_dict = json.loads(params_str)
                                if 'submit_location' in params_dict:
                                    if params_dict['submit_location'] == submit_location:
                                        target_collection = coll
                                        break
                            except:
                                pass
        for i in self.delete_set:
            if target_collection == i[0]:
                delete_targets.append(i)
                break
        if not delete_targets:
            raise Exception('传入的submit_location没有对应的表')
        return delete_targets

    def get_delete_targets(self):
        find_result = self.db['sg_table_relation'].find_one({})
        if find_result:
            target = find_result['target']
        else:
            raise Exception('{}没有"sg_table_relation"这张表'.format(self._project_type))
        return target

    def remove_db_record(self, table_name, query_dict=None, **kwargs):
        """
        根据kwargs删除table_name中查询到的相关记录
        :param table_name: 要查询的表名，如sg_exp
        :param kwargs: 查询条件， 如 diff_id = ObjectId('xxx'), gene_id='ABC'
        :param quey_dict: dict info for querying
        :return:
        """
        conn = self.db[table_name]
        if query_dict:
            kwargs.update(query_dict)
        result = conn.find_one(kwargs)
        if result:
            conn.delete_many(kwargs)
            print('Success to delete records in {} by query {}'.format(table_name, kwargs))
        else:
            print('No record to delete from {} by query {}'.format(table_name, kwargs))

    def remove_table_by_task_id(self, tuple_arg):
        """
        根据task_id删除指定的表，并且删除带有该记录_id的详情表
        :param main_table: 要删除的表的名称，如sg_exp
        :param task_id: task_id对应的值，是删除记录的条件
        :param detail_table: 主表关联的详情表名称如sg_exp_detail，如果有多个详情表，可以是列表,如[sg_xx_1, sg_xx_2]
        :param detail_table_key: 指示是那个字段(如exp_id)对应的值为主表记录的_id, 是删除记录的重要依据。可以是列表，与detail_table一一对应。
        :return:
        """
        main_table, detail_table, detail_table_key = tuple_arg
        all_collections = self.db.collection_names()
        if main_table not in list(all_collections):
            print("Warning: {} was not found in {}".format(main_table, self._project_type))
            return
        conn = self.db[main_table]
        a = conn.find({"task_id": self.task_id, "status": self.status}, {'_id': 1})
        if type(a) == dict:
            a = [a]
        # delete main table record
        remove_num = 0
        for each in a:
            remove_num += 1
            self.remove_db_record(main_table, _id=each['_id'])
            if detail_table:
                if not type(detail_table) == list:
                    detail_table = [detail_table]
                    detail_table_key = [detail_table_key]
                else:
                    if not type(detail_table_key) == list:
                        detail_table_key = [detail_table_key]*len(detail_table)
                for table, table_key in zip(detail_table, detail_table_key):
                    if not table_key:
                        raise Exception('you should specify detail_table_key whose value is main table "_id"')
                    self.remove_db_record(table, query_dict={table_key: each['_id']})
        if remove_num:
            print('Found {} main table(s) in {} by task_id {}.'.format(remove_num, main_table, self.task_id))
            print('And, finished to remove records of {} and {} by task_id {}'.format(main_table, detail_table, self.task_id))

    def run(self):
        self.remove_sg_status()
        target_collections = self.delete_targets
        if not target_collections:
            print('delete_targets为空，结束运行')
            return
        # print('we will delete the following tables:', target_collections)
        with ThreadPoolExecutor(6) as pool:
            pool.map(self.remove_table_by_task_id, target_collections)


if __name__ == '__main__':
    import sys

    if len(sys.argv) < 5:
        exit("Usage: python delete_demo.py <task_id> <project_type> <submit_location> <status>")
    submit_location =sys.argv[3]
    project_type = sys.argv[2]
    task_id = sys.argv[1]
    status = sys.argv[4]
    start_time = time.time()
    del_demo = DeleteRecordsMongo(task_id, project_type, submit_location, status)
    del_demo.run()
    end_time = time.time()
    print("total time: {}".format(end_time - start_time))


