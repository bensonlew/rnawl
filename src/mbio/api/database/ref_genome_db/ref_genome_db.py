# !/usr/bin/python
# -*- coding: utf-8 -*-
from biocluster.api.database.base import Base, report_check
import pymongo
__author__ = 'scp'


class RefGenomeDb(Base):
    def __init__(self, bind_object):
        super(RefGenomeDb, self).__init__(bind_object)
        self._project_type = 'ref_rna_v2'

    def insert_document(self, table_name, content_dict_list, client, tag_dict=None):
        table_id = None
        if client == "client03":
            conn = self.db[table_name]
        else:
            tsg_client = pymongo.MongoClient(
                "mongodb://rna:y6a3t5n1w0y7@10.8.0.23/sanger_ref_rna_v2?authMechanism=SCRAM-SHA-1")
            db = tsg_client["sanger_ref_rna_v2"]
            conn = db[table_name]
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

    def remove_db_record(self, table_name, query_dict=None, **kwargs):
        """
        Delete related records queried in table_name according to kwargs
        :param table_name: table name, such as sg_exp
        :param kwargs: query conditions, such as diff_id=ObjectId('5c04df94edcb2573e7b68169')
        :param quey_dict: dict info for querying
        :return:
        """
        client = self.bind_object.sheet.client
        if client == "client03":
            conn = self.db[table_name]
        else:
            tsg_client = pymongo.MongoClient(
                "mongodb://rna:y6a3t5n1w0y7@10.8.0.23/sanger_ref_rna_v2?authMechanism=SCRAM-SHA-1")
            db = tsg_client["sanger_ref_rna_v2"]
            conn = db[table_name]
        if query_dict:
            kwargs.update(query_dict)
        result = conn.find_one(kwargs)
        if result:
            conn.delete_many(kwargs)
        else:
            print('no record to delete')