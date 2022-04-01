# !/usr/bin/python
# -*- coding: utf-8 -*-
from biocluster.api.database.base import Base, report_check
import types
from bson.objectid import ObjectId
from collections import OrderedDict
__author__ = 'gdq'


class ApiBase(Base):
    def __init__(self, bind_object):
        super(ApiBase, self).__init__(bind_object)
        self._project_type = 'prok_rna'

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
                print("Failed to insert records into table {} as: {}".format(table_name, e))
            else:
                print("Failed to insert record into table {} as: {}".format(table_name, e))
        else:
            if record_num >= 2:
                print("Success to insert records into table {}".format(table_name))
            else:
                print("Success to insert record into table {}".format(table_name))

            return table_id

    def update_snpid_2objected(self, table_name):
        conn = self.db[table_name]
        conn.update()

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
            print('No record to delete')

    def update_db_record(self, table_name, record_id=None, query_dict=None, insert_dict=None, **kwargs):
        if record_id is not None:
            if isinstance(record_id, types.StringTypes):
                record_id = ObjectId(record_id)
            elif isinstance(record_id, ObjectId):
               record_id = record_id
            else:
                raise Exception("main_id参数必须为字符串或者ObjectId类型!")
        conn = self.db[table_name]
        if query_dict:
            if record_id is not None:
                query_dict.update({"_id": record_id})
        else:
            if record_id is not None:
                query_dict = {"_id": record_id}
            else:
                raise Exception("Please provide query dict")
        if insert_dict:
            kwargs.update(insert_dict)
        conn.update(query_dict, {"$set": kwargs}, upsert=True)

    def update_db_record_2(self, table_name, record_id, query_dict=None, **kwargs):
        if isinstance(record_id, types.StringTypes):
            record_id = ObjectId(record_id)
        elif isinstance(record_id, ObjectId):
           record_id = record_id
        else:
            raise Exception("main_id参数必须为字符串或者ObjectId类型!")
        conn = self.db[table_name]
        if query_dict:
            kwargs.update(query_dict)
        conn.update({"_id": record_id}, {"$set": kwargs}, upsert=True)

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

    def update_sgtask_record(self, task_id=None, table_name="sg_task", **kwargs):
        """
        根据记录的唯一标识tsk_id找到sg_task的相关记录，并用kwargs更新记录
        :param table_name:
        :param task_id:
        :param kwargs: kwargs to add
        :return:
        """
        if task_id is None:
            if self.bind_object:
                task_id=self.bind_object.sheet.id
            else:
                raise Exception('No valid task_id')
        conn = self.db[table_name]
        conn.update({"task_id": task_id}, {"$set": kwargs}, upsert=True)

    def get_project_sn(self, task_id=None):
        if task_id is None:
            if self.bind_object:
                task_id = self.bind_object.sheet.id
            else:
                raise Exception('No valid task_id')
        conn = self.db["sg_task"]
        result = conn.find_one({"task_id": task_id})
        return result['project_sn']

    def update_sample_order(self, exp_id, group_dict, ):
        """
        暂时无用
        :param exp_id:
        :param group_dict:
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
        根据task_id删除指定的表，并且删除带有该记录_id的详情表
        :param main_table: 要删除的表的名称，如sg_exp
        :param task_id: task_id对应的值，是删除记录的条件
        :param detail_table: 主表关联的详情表名称如sg_exp_detail，如果有多个详情表，可以是列表,如[sg_xx_1, sg_xx_2]
        :param detail_table_key: 指示是那个字段(如exp_id)对应的值为主表记录的_id, 是删除记录的重要依据。可以是列表，与detail_table一一对应。
        :return:
        """
        conn = self.db[main_table]
        a = conn.find({"task_id": task_id }, {'_id': 1})
        if type(a) == dict():
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
        print('Have removed records of {} and {} by {}'.format(main_table, detail_table, task_id))

    def remove_table_by_main_record(self, main_table, detail_table=None, detail_table_key=None, **kwargs):
        """
        :根据maintable 属性 删除指定的表，并且删除带有该记录_id的详情表
        :param main_table: 要删除的表的名称，如sg_exp
        :param task_id: task_id对应的值，是删除记录的条件
        :param detail_table: 主表关联的详情表名称如sg_exp_detail，如果有多个详情表，可以是列表,如[sg_xx_1, sg_xx_2]
        :param detail_table_key: 指示是那个字段(如exp_id)对应的值为主表记录的_id, 是删除记录的重要依据。可以是列表，与detail_table一一对应。
        :return:
        :param kwargs 其它属性表示匹配需要删除主表的属性
        """
        conn = self.db[main_table]
        a = conn.find(kwargs, {'_id': 1})
        if type(a) == dict():
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
        print('Have removed records of {} and {} by {}'.format(main_table, detail_table, kwargs))

    def get_table_by_main_record(self, main_table, **kwargs):
        """
        :根据maintable 属性获取表格记录
        :param main_table: 要删除的表的名称，如sg_exp
        :param kwargs 其它属性表示匹配需要删除主表的属性
        """
        conn = self.db[main_table]
        a = conn.find_one(kwargs, {'_id': 1})
        return a

    def order_row_dict_list(self, row_dict_list, order_list):
        """
        对pandas的to_dict函数的结果进行排序
        :param row_dict_list:
        :param order_list: 排序依据，还可以借此过滤出想要的列信息。
        :return:
        """
        ordered_dict_list = list()
        for each in row_dict_list:
            tmp_dict = OrderedDict()
            for sample in order_list:
                if sample in each:
                    tmp_dict[sample] = each[sample]
                else:
                    print("Sample named {} is not in order list".format(sample))
            ordered_dict_list.append(tmp_dict)
        return ordered_dict_list
