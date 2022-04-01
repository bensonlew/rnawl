# -*- coding: utf-8 -*-
# __author__ = 'guoquan'  
from mainapp.config.db import Config
import web
import re


class BaseModel(object):

    def __init__(self, table=None):
        self._db = Config().get_db()
        self._table = table

    @property
    def db(self):
        return self._db

    @property
    def table(self):
        if not self._table:
            self._table = self._get_table_name()
        return self._table

    def _get_table_name(self):
        class_name = self.__class__.__name__
        return re.sub(r"Model", "", class_name).lower()

    def get_by_id(self, keyid, what="*"):
        where_str = "id = %s" % keyid
        data = self.db.select(self.table, what=what, where=where_str)
        if len(data) > 0:
            return data[0]
        else:
            return None


class BaseListModel(BaseModel):

    def __init__(self, table=None):
        super(BaseListModel, self).__init__(table)
        self._pagesize = 100
        self._table = table
        self._total = None
        self._where = None
        self._what = "*"
        self._order = "id DESC"
        self._current_page = 1

    @property
    def current_page(self):
        return self._current_page

    def set_conditions(self, where=None, what="*", order="id DESC"):
        self._where = where
        self._what = what
        self._order = order

    @property
    def total(self):
        if self._total is None:
            data = self.db.select(self.table, what="count(id) as total", where=self._where)
            self._total = data[0].total
        return self._total

    @property
    def total_page(self):
        page = self.total / self.pagesize
        if self.total > page * self.pagesize:
            page += 1
        return page

    @property
    def pagesize(self):
        return self._pagesize

    @pagesize.setter
    def pagesize(self, size):
        if size < 1:
            raise Exception("size must more then 1")
        self._pagesize = size

    def get_data(self, page=1):
        if page < 1:
            page = 1
        if page > self.total_page:
            page = self.total_page
        self._current_page = page
        offset = (page - 1) * self._pagesize
        data = self.db.select(self.table, what=self._what, where=self._where,
                              order=self._order, limit=self.pagesize, offset=offset)
        return data



