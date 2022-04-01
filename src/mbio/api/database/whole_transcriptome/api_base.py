# -*- coding: utf-8 -*-
# __author__ = 'gudeqing'

import types
from collections import OrderedDict

from biocluster.api.database.base import Base
from bson.objectid import ObjectId


class ApiBase(Base):
    def __init__(self, bind_object):
        super(ApiBase, self).__init__(bind_object)
        self._project_type = 'whole_transcriptome'

    def create_db_table(self, table_name, content_dict_list, tag_dict=None):
        '''
        Create main/detail table in database system.
        :param table_name: table name
        :param content_dict_list: list with dict as elements
        :param tag_dict: a dict to be added into each record in content_dict_list.
        :return: None or main table id
        '''
        table_id = None
        conn = self.db[table_name]
        if tag_dict:
            for row_dict in content_dict_list:
                row_dict.update(tag_dict)
        record_num = len(content_dict_list)
        try:
            if record_num > 5000:
                for i in range(0, record_num, 3000):
                    tmp_list = content_dict_list[i: i + 3000]
                    conn.insert_many(tmp_list)
            else:
                if record_num >= 2:
                    conn.insert_many(content_dict_list)
                else:
                    table_id = conn.insert_one(content_dict_list[0]).inserted_id
        except Exception as e:
            if record_num >= 2:
                self.bind_object.logger.warn('fail to insert records into table {} -> ({})'.format(table_name, e))
            else:
                self.bind_object.logger.warn('fail to insert record into table {} -> ({})'.format(table_name, e))
        else:
            if record_num >= 2:
                self.bind_object.logger.info('succeed in inserting records into table {}'.format(table_name))
            else:
                self.bind_object.logger.info('succeed in inserting record into table {}'.format(table_name))
            return table_id

    def remove_db_record(self, table_name, query_dict=None, **kwargs):
        '''
        Delete related records queried in table_name according to kwargs
        :param table_name: the name of the table to be queried, e.g. sg_exp
        :param kwargs: query conditions, e.g. diff_id=ObjectId('5d7f1f65ffec602b062b746d'), gene_id='ENSG00000109062'
        :param quey_dict: dict info for querying
        :return:
        '''
        conn = self.db[table_name]
        if query_dict:
            kwargs.update(query_dict)
        result = conn.find_one(kwargs)
        if result:
            conn.delete_many(kwargs)
        else:
            self.bind_object.logger.warn('fail to find {} by {}'.format(table_name, kwargs))

    def update_db_record(self, table_name, record_id=None, query_dict=None, insert_dict=None, **kwargs):
        if record_id is not None:
            assert isinstance(record_id, types.StringTypes) or isinstance(record_id, ObjectId)
            record_id = ObjectId(record_id)
        conn = self.db[table_name]
        if query_dict:
            if record_id is not None:
                query_dict.update({'_id': record_id})
        else:
            if record_id is not None:
                query_dict = {'_id': record_id}
        if insert_dict:
            kwargs.update(insert_dict)
        conn.update(query_dict, {'$set': kwargs}, upsert=True)

    def update_db_record_by_dict(self, table_name, record_dict, **kwargs):
        """
        根据记录的标识dict找到table_name的相关记录，并用kwargs更新记录
        :param table_name:
        :param record_dict:
        :param kwargs: kwargs to add
        :return:
        """
        conn = self.db[table_name]
        conn.update(record_dict, {"$set": kwargs}, upsert=True)

    def update_sgtask_record(self, task_id=None, table_name='sg_task', **kwargs):
        if task_id is None:
            task_id = self.bind_object.sheet.id
        conn = self.db[table_name]
        conn.update({'task_id': task_id}, {'$set': kwargs}, upsert=True)

    def get_project_sn(self, task_id=None):
        if task_id is None:
            task_id = self.bind_object.sheet.id
        conn = self.db['sg_task']
        result = conn.find_one({'task_id': task_id})
        return result['project_sn']

    def update_record_by_task_id(self, table_name, task_id, **kwargs):
        conn = self.db[table_name]
        conn.update({'task_id': task_id}, {'$set': kwargs}, upsert=True)

    def remove_table_by_task_id(self, main_table, task_id, detail_table=None, detail_table_key=None):
        conn = self.db[main_table]
        main_records = conn.find({'task_id': task_id}, {'_id': 1})
        for main_record in main_records:
            self.remove_db_record(main_table, _id=main_record['_id'])
            if detail_table:
                if type(detail_table) == str:
                    detail_table = [detail_table]
                    detail_table_key = [detail_table_key]
                else:
                    if type(detail_table_key) == str:
                        detail_table_key = [detail_table_key] * len(detail_table)
                for table, table_key in zip(detail_table, detail_table_key):
                    self.remove_db_record(table, query_dict={table_key: main_record['_id']})
        else:
            self.bind_object.logger.info('succeed in removing tables in {} and {} by {}'.format(main_table, detail_table, task_id))

    def remove_table_by_main_record(self, main_table, detail_table=None, detail_table_key=None, **kwargs):
        conn = self.db[main_table]
        main_records = conn.find(kwargs, {'_id': 1})
        if type(main_records) == dict:
            main_records = [main_records]
        for main_record in main_records:
            self.remove_db_record(main_table, _id=main_record['_id'])
            if detail_table:
                if type(detail_table) == str:
                    detail_table = [detail_table]
                    detail_table_key = [detail_table_key]
                else:
                    if type(detail_table_key) == str:
                        detail_table_key = [detail_table_key] * len(detail_table)
                for table, table_key in zip(detail_table, detail_table_key):
                    self.remove_db_record(table, query_dict={table_key: main_record['_id']})
        else:
            self.bind_object.logger.info('succeed in removing tables in {} and {} by {}'.format(main_table, detail_table, kwargs))

    def get_table_by_main_record(self, main_table, **kwargs):
        conn = self.db[main_table]
        document = conn.find_one(kwargs, {'_id': 1})
        return document

    def order_row_dict_list(self, row_dict_list, order_list):
        ordered_dict_list = list()
        for each in row_dict_list:
            tmp_dict = OrderedDict()
            for sample in order_list:
                if sample in each:
                    tmp_dict[sample] = each[sample]
            ordered_dict_list.append(tmp_dict)
        return ordered_dict_list
