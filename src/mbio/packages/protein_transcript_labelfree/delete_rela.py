# -*- coding: utf-8 -*-
# __author__ = 'gdq'
import time
from biocluster.config import Config
from concurrent.futures import ThreadPoolExecutor


class DeleteRelaMongo(object):
    def __init__(self, task_id, project_type, delete_targets=None):
        self.db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        self.task_id = task_id
        self._project_type = project_type
        if not delete_targets:
            self.delete_targets = self.get_delete_targets()
        else:
            self.delete_targets = delete_targets

    def get_delete_targets(self):
        find_result = self.db['sg_table_relation_for_rel'].find_one({})
        if find_result:
            target = find_result['target']
        else:
            raise Exception('{}没有"sg_table_relation_for_rel"这张表'.format(self._project_type))
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
        a = conn.find({"task_id": self.task_id}, {'_id': 1})
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
        target_collections = self.delete_targets
        if not target_collections:
            print('delete_targets为空，结束运行')
            return
        # print('we will delete the following tables:', target_collections)
        with ThreadPoolExecutor(6) as pool:
            pool.map(self.remove_table_by_task_id, target_collections)
        self.deal_sg_status()

    def deal_sg_status(self):
        conn = self.db['sg_status']
        # print('111')
        a = conn.find({"task_id": self.task_id}, {'_id': 1, 'submit_location':1})
        # print('222')
        # print(list(a))
        if type(a) == dict:
            a = [a]
        for each in a:
            print('111')
            if u'submit_location' in each and each['submit_location'] in [ 'p2grelationship', 'proteingene', 'p2gdiff', 'relasetcluster', 'relasetcorr', 'prcorr', 'delete_relaset', 'delete_rela', 'prgo', 'prgorich', 'richcluster', 'prkegg', 'prkeggrich', 'keggrichcluster']:
                print(each['submit_location'])
            #self.remove_db_record('sg_status', _id=each['_id'])
                conn.update({'_id': each['_id']}, {'$set': {'status': 'deleted'}})

if __name__ == '__main__':
    import sys

    print("Usage: python delete_rela.py <task_id> <project_type>")
    project_type = sys.argv[2]
    task_id = sys.argv[1]
    start_time = time.time()
    del_demo = DeleteRelaMongo(task_id, project_type, delete_targets=None)
    del_demo.run()
    end_time = time.time()
    print("total time: {}".format(end_time - start_time))


