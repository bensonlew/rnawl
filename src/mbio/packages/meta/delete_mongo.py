# -*- coding: utf-8 -*-
# __author__ = 'gudeqing'

import logging
import sys
import time

from biocluster.config import Config
from concurrent.futures import ThreadPoolExecutor

# logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
#                     datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)


class DeleteDemoMongo(object):
    def __init__(self, task_id, project_type, delete_targets=None, db_version=None):
        if db_version:
            self.db = Config().get_mongo_client(mtype=project_type, db_version=db_version)[Config().get_mongo_dbname(mtype=project_type, db_version=db_version)]
            self.db_static = Config().get_mongo_client(mtype=project_type, db_version=db_version, dydb_forbid=True)[
                Config().get_mongo_dbname(mtype=project_type, db_version=db_version, dydb_forbid=True)]
        else:
            self.db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
            self.db_static = Config().get_mongo_client(mtype=project_type, dydb_forbid=True)[Config().get_mongo_dbname(project_type, dydb_forbid=True)]
        self.task_id = task_id
        self._project_type = project_type
        if not delete_targets:
            self.delete_targets = self.get_delete_targets()
        else:
            self.delete_targets = delete_targets

    def get_delete_targets(self):
        """
        if self._project_type == 'whole_transcriptome':
            document = self.db_static['table_relation'].find_one({})
            if document and 'target' in document:
                return self.db_static['table_relation'].find_one({})['target']
        """
        project_type = ['metagenomic', "metagbin"]
        if self._project_type in project_type:
            documents = self.db_static['table_relation'].find({})
            target = []
            for document in documents:
                every_target = []
                if str(document["is_detail"]) == "y":
                    every_target.append(document["analysis"])
                    every_target.append(document["table"])
                    every_target.append(document["asso_id"])
                if every_target:
                    target.append(every_target)
        else:
            find_result = self.db_static['sg_table_relation'].find_one({})
            if find_result:
                target = find_result['target']
            else:
                raise Exception('can not find sg_table_relation with project type ({})'.format(self._project_type))

        target_dict = {}
        for i in target:
            main_table = i[0]
            detail_table = i[1]
            relation_id = i[2]
            if main_table not in target_dict.keys():
                target_dict[main_table] = {}
                if relation_id not in target_dict[main_table].keys():
                    target_dict[main_table][relation_id] = []
                    if detail_table not in target_dict[main_table][relation_id]:
                        target_dict[main_table][relation_id].append(detail_table)
            else:
                if relation_id not in target_dict[main_table].keys():
                    target_dict[main_table][relation_id] = []
                    if detail_table not in target_dict[main_table][relation_id]:
                        target_dict[main_table][relation_id].append(detail_table)
                else:
                    if detail_table not in target_dict[main_table][relation_id]:
                        target_dict[main_table][relation_id].append(detail_table)

        new_target = []
        for i in target_dict.keys():
            for j in target_dict[i]:
                tmp_target = []
                tmp_target.append(i)
                tmp_target.append(target_dict[i][j])
                tmp_target.append(j)
                new_target.append(tmp_target)
        return new_target

    def run(self):
        target_collections = self.delete_targets
        if not target_collections:
            logging.info('targets for deleting is empty, abord')
            return
        with ThreadPoolExecutor(8) as pool:
            pool.map(self.remove_table_by_task_id, target_collections)
        self.deal_sg_status()

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
            logging.warn('{} was not found with project type ({})'.format(main_table, self._project_type))
            return
        if main_table in ["raw_stat_detail", "qc_stat_detail"]:  # 这里是一个特殊情况，bac_assem流程中的graphic_detail详情表关联的不是主表，而是这两张详情表
            cursor = []
            tmp_main_table = self.db["sample_info"].find({'task_id': self.task_id}, {'_id': 1})  # 通过task_id在sample_info这张主表里查找到关联的所有这两张详情表，作为graphic_detail的主表
            for table in tmp_main_table:
                relation_id = table["_id"]
                tmp_cursor = self.db[main_table].find({'info_id': relation_id}, {'_id': 1})
                for a_cursor in tmp_cursor:
                    cursor.append(a_cursor)
        else:
            collection = self.db[main_table]
            cursor = collection.find({'task_id': self.task_id}, {'_id': 1})
        if type(cursor) == dict:
            cursor = [cursor]
        remove_num = 0
        for main_document in cursor:
            remove_num += 1
            self.remove_db_record(main_table, _id=main_document['_id'])
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
                    self.remove_db_record(table, query_dict={table_key: main_document['_id']})
        if remove_num:
            logging.info(
                'found {} main table(s) in {} by task_id ({}), removed records of {}'.format(remove_num, main_table,
                                                                                             self.task_id,
                                                                                             detail_table))

    def remove_db_record(self, table_name, query_dict=None, **kwargs):
        """
        根据kwargs删除table_name中查询到的相关记录
        :param table_name: 要查询的表名，如sg_exp
        :param kwargs: 查询条件， 如 diff_id = ObjectId('xxx'), gene_id='ABC'
        :param quey_dict: dict info for querying
        :return:
        """
        collection = self.db[table_name]
        if query_dict:
            kwargs.update(query_dict)
        document = collection.find_one(kwargs)
        if document:
            collection.delete_many(kwargs)
            logging.info('succeed in deleting records in {} by query {}'.format(table_name, kwargs))
        else:
            logging.info('find no record to delete from {} by query {}'.format(table_name, kwargs))

    def deal_sg_status(self):
        collection = self.db['sg_status']
        cursor = collection.find({"task_id": self.task_id}, {'_id': 1})
        if type(cursor) == dict:
            cursor = [cursor]
        for document in cursor:
            collection.update({'_id': document['_id']}, {'$set': {'status': 'deleted'}})


if __name__ == '__main__':
    logging.info('Usage: python delete_demo.py <task_id> <project_type> <db_version>')
    task_id = sys.argv[1]
    project_type = sys.argv[2]
    if len(sys.argv) >= 4:
        db_version = sys.argv[3]
    start_time = time.time()
    inst = DeleteDemoMongo(task_id, project_type, delete_targets=None, db_version=db_version)
    inst.run()
    end_time = time.time()
    logging.info('Elapsed time: {}'.format(end_time - start_time))
