# -*- coding: utf-8 -*-
# __author__ = ''
from .update_status import UpdateStatus


class TupdateStatus(UpdateStatus):
    def __init__(self, data):
        super(TupdateStatus, self).__init__(data)
        self._client = "client03"
        # self._key = "hM4uZcGs9d"
        # self._url = "http://api.tsg.com/task/add_file"
        # self._url = "http://api.tsanger.com/task/add_file"
        self._project_type = 'metabolome'
        self._env_name = "offline"
        self._key = "458b97de4c0bb5bf416c8cea208309ed"
        self._url = "http://apicenter.nsg.com/index/in"
