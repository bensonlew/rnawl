# -*- coding: utf-8 -*-
# __author__ = 'guoquan'  

from mainapp.config.db import Config
import web
import datetime


class User(object):

    def __init__(self, username=None):
        self._db = Config().get_db()
        self._username = username
        self._data = None
        self.table = "user"

    @property
    def data(self):
        if not self._data:
            self._data = self._get_info()
        return self._data

    def _get_info(self):
        where_dict = dict(user=self._username)
        data = self._db.select(self.table, where=web.db.sqlwhere(where_dict))
        if len(data) > 0:
            return data[0]
        else:
            return False

    def check_pass(self, password):
        if not self.data:
            return False
        if self.data.lock:
            return False
        if self.data.password == password:
            return True
        return False

    def is_admin(self):
        if not self.data:
            return False
        if self.data.lock:
            return False
        if self.data.supper == 1:
            return True
        return False

    def login(self):
        self._db.update(self.table, where=dict(user=self._username), last_login=datetime.datetime.now(),
                        last_ip=web.ctx.ip)
