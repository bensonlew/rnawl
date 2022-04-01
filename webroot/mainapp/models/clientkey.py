# -*- coding: utf-8 -*-
# __author__ = 'HD'
from mainapp.config.db import Config
import web


class ClientKey(object):
    def __init__(self, client):
        self.config = Config()
        self.client = client
        self.key = None
        self.ipmask = None
        self.timelimit = None
        self.n_key = None
        self._select()
        self._select_new_key()

    def _select(self):
        if self.client not in self.config.CON_FILE_DICT.keys():
            raise web.badrequest("client not found!")
        else:
            self.key = self.config.CON_FILE_DICT[self.client]["key"]
            self.ipmask = self.config.CON_FILE_DICT[self.client]["ipmask"]
            self.timelimit = self.config.CON_FILE_DICT[self.client]["timelimit"]

    def _select_new_key(self):
        """
        add by hd @20200821 获取mysql中的n_client对应的key，该key是前端验证的新key，以前的client中的key也要保留作用
        :return:
        """
        if self.client not in ["client01", "client03"]:
            return
        if "n_" + self.client not in self.config.CON_FILE_DICT.keys():
            raise web.badrequest("client not found!")
        else:
            self.n_key = self.config.CON_FILE_DICT["n_" + self.client]["key"]
