# -*- coding: utf-8 -*-
# __author__ = 'fyt'
import time
from biocluster.api.database.base import Base
import gevent.pool
import gevent.monkey
import json
from bson import CodecOptions, SON

class DeleteMongo(Base):
    def __init__(self, project_type='prok_rna'):
        super(DeleteMongo, self).__init__()
        self._project_type = project_type
        self.db_name = self.get_db_name()

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
        else:
            print('No record to delete from {} by {}'.format(table_name, kwargs))

    def get_db_name(self):
        target_collections = []
        collections = self.db.collection_names()
        for coll in collections:
            # if coll.endswith('_detail') and coll != 'sg_t2g_detail':
            if coll.endswith('_detail'):
                continue
            find_result = self.db[coll].find_one({'params': {'$exists': True}})
            if find_result:
                target_collections.append(coll)  # 过滤出含有params字段的主表
        return target_collections

    def judge_db_null(self):
        params_null = False
        collections = self.db.collection_names()
        for coll in collections:
            # if coll.endswith('_detail') and coll != 'sg_t2g_detail':
            if coll.endswith('_detail'):
                continue
            find_result = self.db[coll].find_one({'params': {'$exists': True}})
            if find_result:
                params_null = True  # 过滤出含有params字段的主表
                break
        return params_null

    def delete_mongo(self,tuple_arg):
        main_table, detail_table, detail_table_key = tuple_arg
        if main_table in self.db_name:
            # print(main_table)
            conn = self.db[main_table]
            a = conn.find({'params': {'$exists': True}})
            if type(a) == dict:
                a = [a]
            # delete main table record
            for each in a:
                if each['params'] is None:
                    self.remove_db_record(main_table, _id=each['_id'])
                    if detail_table:
                        if not type(detail_table) == list:
                            detail_table = [detail_table]
                            detail_table_key = [detail_table_key]
                        else:
                            if not type(detail_table_key) == list:
                                detail_table_key = [detail_table_key] * len(detail_table)
                        for table, table_key in zip(detail_table, detail_table_key):
                            if not table_key:
                                raise Exception('you should specify detail_table_key whose value is main table "_id"')
                            self.remove_db_record(table, query_dict={table_key: each['_id']})
                    print('{} We removed records of {} and {},the objectid is {}\n'.format(time.strftime("%Y%m%d_%H%M%S"), main_table, detail_table, each['_id']))

    def get_delete_targets(self):
        target = []
        try:
            find_result = self.db['sg_table_relation'].find_one({})
            if find_result:
                target = find_result['target']
        except:
            print('{}没有"sg_table_relation"这张表'.format(self._project_type))
            target = []
        return target

    def run(self):
        delete_target=self.get_delete_targets()
        # print(delete_target[0])
        pool = gevent.pool.Pool(100)
        pool.map(self.delete_mongo, delete_target)

if __name__ == '__main__':
    projects = ['prok_rna', 'denovo_rna_v2', 'itraq_tmt', 'labelfree']
    delete = DeleteMongo()
    del_cmd=1
    while del_cmd:
        pm=time.strftime('%H')
        if pm != '12':
            time.sleep(3600)
        else:
            print("\n老铁，现在是{}，开始干活\n".format(time.strftime("%Y%m%d_%H%M%S")))
            time_start=time.time()
            time_end = 0
            for project in projects:
                delete = DeleteMongo(project)
                if not delete.judge_db_null():
                    print('{} 今天在{}项目中没有找到可以清除的表'.format(time.strftime("%Y%m%d_%H%M%S"), project))
                    time_end = time.time()
                else:
                    print('\n开始清除{}'.format(project))
                    delete.run()
                    print('\n已经清除完{}'.format(project))
                    time_end=time.time()
            time.sleep(3600 - time_end + time_start)
