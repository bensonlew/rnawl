# !/usr/bin/python
# -*- coding: utf-8 -*-
import re
from collections import defaultdict
import sqlite3
import datetime
from biocluster.api.database.base import Base, report_check
import unittest

class DeleteDemo(Base):
    def __init__(self, bind_object):
        super(DeleteDemo, self).__init__(bind_object)
        self._project_type = 'project'

    def update_main_table(self, collection_name, main_id, status):
        """
        更新主表字段
        :param collection_name:
        :param data:
        :return: 插入记录的id
        """
        conn = self.db[collection_name]
        conn.update({'_id': main_id}, {"$set": {'del_mongo_status': status}}, upsert=True)