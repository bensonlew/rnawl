# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from ..meta.update_status import UpdateStatus


class TupdateStatus(UpdateStatus):

    def __init__(self, data):
        super(TupdateStatus, self).__init__(data)
        self._client = "client03"
        # self._key = "hM4uZcGs9d"
        # self._url = "http://www.tsanger.com/api/add_file"
        self.mongodb = self._mongo_client[self.config.MONGODB + '_rna']
        self._env_name = "offline"
        self._key = "458b97de4c0bb5bf416c8cea208309ed"
        self._url = "http://apicenter.nsg.com/index/in"
