# -*- coding: utf-8 -*-
# __author__ = 'gudeqing'

from __future__ import print_function
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
import types
from collections import OrderedDict

class ApiBase(Base):
    def __init__(self, bind_object):
        super(ApiBase, self).__init__(bind_object)
        self._project_type = 'lnc_rna'

    def create_db_table(self, table_name, content_dict_list, tag_dict=None):
        """
        Create main/detail table in database system.
        :param table_name: table name
        :param content_dict_list: list with dict as elements
        :param tag_dict: a dict to be added into each record in content_dict_list.
        :return: None or main table id
        """
        table_id = None
        conn = self.db[table_name]
        if tag_dict:
            for row_dict in content_dict_list:
                row_dict.update(tag_dict)
        record_num = len(content_dict_list)
        try:
            # split data and dump them to db separately
            if record_num > 5000:
                for i in range(0, record_num, 3000):
                    tmp_list = content_dict_list[i: i+3000]
                    conn.insert_many(tmp_list)
            else:
                if record_num >= 2:
                    conn.insert_many(content_dict_list)
                else:
                    table_id = conn.insert_one(content_dict_list[0]).inserted_id
        except Exception as e:
            if record_num >= 2:
                print("fail to insert records into table {} as: {}".format(table_name, e))
            else:
                print("fail to insert record into table {} as: {}".format(table_name, e))
        else:
            if record_num >= 2:
                print("succeed in inserting records into table {}".format(table_name))
            else:
                print("succeed in inserting record into table {}".format(table_name))
            return table_id

    def remove_db_record(self, table_name, query_dict=None, **kwargs):
        """
        Delete related records queried in table_name according to kwargs
        :param table_name: table name, such as sg_exp
        :param kwargs: query conditions, such as diff_id=ObjectId('5c04df94edcb2573e7b68169')
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
            print('no record to delete')

    def update_db_record(self, table_name, record_id=None, query_dict=None, insert_dict=None, **kwargs):
        if record_id is not None:
            if isinstance(record_id, types.StringTypes):
                record_id = ObjectId(record_id)
            elif isinstance(record_id, ObjectId):
               record_id = record_id
            else:
                raise Exception('main_id must be a string or ObjectId type')
        conn = self.db[table_name]
        if query_dict:
            if record_id is not None:
                query_dict.update({"_id": record_id})
        else:
            if record_id is not None:
                query_dict = {"_id": record_id}
            else:
                raise Exception("please provide record id or query dict")
        if insert_dict:
            kwargs.update(insert_dict)
        conn.update(query_dict, {"$set": kwargs}, upsert=True)

    def update_sgtask_record(self, task_id=None, table_name="sg_task", **kwargs):
        """
        Find sg_task related record based on task_id and update the record with kwargs
        :param table_name: table name
        :param task_id: task id
        :param kwargs: kwargs to add
        :return:
        """
        if task_id is None:
            if self.bind_object:
                task_id=self.bind_object.sheet.id
            else:
                raise Exception('no valid task_id')
        conn = self.db[table_name]
        conn.update({"task_id": task_id}, {"$set": kwargs}, upsert=True)

    def get_project_sn(self, task_id=None):
        if task_id is None:
            if self.bind_object:
                task_id = self.bind_object.sheet.id
            else:
                raise Exception('no valid task_id')
        conn = self.db["sg_task"]
        result = conn.find_one({"task_id": task_id})
        ret = None
        if 'project_sn' in result:
            ret = result['project_sn']
        else:
            print('no valid project_sn')
        return ret

    def update_sample_order(self, exp_id, group_dict, ):
        """
        Deprecated
        :param exp_id: exp id
        :param group_dict: group dict
        :return:
        """
        sample_order = list()
        group_order = list()
        for k in group_dict:
            sample_order += group_dict[k]
            group_order.append(k)
        self.update_db_record('sg_exp', exp_id, sample_order=sample_order, group_order=group_order)

    def update_record_by_task_id(self, table_name, task_id, **kwargs):
        conn = self.db[table_name]
        conn.update({"task_id": task_id}, {"$set": kwargs}, upsert=True)

    def remove_table_by_task_id(self, main_table, task_id, detail_table=None, detail_table_key=None):
        """
        Delete specified main table based on task_id and detail table with related _id
        :param main_table: table name, such as sg_exp
        :param task_id: task id
        :param detail_table: both str and list are supported, such as sg_exp_detail or [sg_diff_detail, sg_diff_scatter]
        :param detail_table_key: key(s) used to associate main table with detail table, must be correspond to detail_table
        :return:
        """
        conn = self.db[main_table]
        a = conn.find({"task_id": task_id }, {'_id': 1})
        if type(a) == dict:
            a = [a]
        # delete main table record
        for each in a:
            self.remove_db_record(main_table, _id=each['_id'])
            if detail_table:
                if type(detail_table) == str:
                    detail_table = [detail_table]
                    detail_table_key = [detail_table_key]
                else:
                    if type(detail_table_key) == str:
                        detail_table_key = [detail_table_key]*len(detail_table)
                for table, table_key in zip(detail_table, detail_table_key):
                    if not table_key:
                        raise Exception('you should specify detail_table_key whose value is main table "_id"')
                    self.remove_db_record(table, query_dict={table_key: each['_id']})
        print('have removed records of {} and {} by {}'.format(main_table, detail_table, task_id))

    def remove_table_by_main_record(self, main_table, detail_table=None, detail_table_key=None, **kwargs):
        """
        Delete specified main table based on task_id and detail table with related _id
        :param main_table: table name, such as sg_exp
        :param detail_table: both str and list are supported, such as sg_exp_detail or [sg_diff_detail, sg_diff_scatter]
        :param detail_table_key: key(s) used to associate main table with detail table, must be correspond to detail_table
        :param kwargs: property of the main table that need to be deleted
        :return:
        """
        conn = self.db[main_table]
        a = conn.find(kwargs, {'_id': 1})
        if type(a) == dict:
            a = [a]
        # delete main table record
        for each in a:
            self.remove_db_record(main_table, _id=each['_id'])
            if detail_table:
                if type(detail_table) == str:
                    detail_table = [detail_table]
                    detail_table_key = [detail_table_key]
                else:
                    if type(detail_table_key) == str:
                        detail_table_key = [detail_table_key]*len(detail_table)
                for table, table_key in zip(detail_table, detail_table_key):
                    if not table_key:
                        raise Exception('you should specify detail_table_key whose value is main table "_id"')
                    self.remove_db_record(table, query_dict={table_key: each['_id']})
        print('have removed records of {} and {} by {}'.format(main_table, detail_table, kwargs))

    def get_table_by_main_record(self, main_table, **kwargs):
        """
        Get the table records based on the property of the main table
        :param main_table: table name, such as sg_exp
        :param kwargs: property of the main table that need to be obtained
        :return: main table
        """
        conn = self.db[main_table]
        return conn.find_one(kwargs, {'_id': 1})

    def get_dict_by_main_record(self, main_table, **kwargs):
        """
        Get the table records based on the property of the main table
        :param main_table: table name, such as sg_exp
        :param kwargs: property of the main table that need to be obtained
        :return: main table
        """
        conn = self.db[main_table]
        return conn.find_one(kwargs)

    def order_row_dict_list(self, row_dict_list, order_list):
        """
        Sort the results from to_dict function of pandas
        :param row_dict_list: list containing dict(s)
        :param order_list: base for sorting
        :return:
        """
        ordered_dict_list = list()
        for each in row_dict_list:
            tmp_dict = OrderedDict()
            for sample in order_list:
                if sample in each:
                    tmp_dict[sample] = each[sample]
                else:
                    print("sample named {} is not in order list".format(sample))
            ordered_dict_list.append(tmp_dict)
        return ordered_dict_list
