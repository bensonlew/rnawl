# -*- coding: utf-8 -*-
# __author__ = 'HD'
from biocluster.config import Config
from .update_status import UpdateStatus


class TupdateStatus(UpdateStatus):

    def __init__(self, data):
        super(TupdateStatus, self).__init__(data)
        # self._client = Config().rcf.get("toollab", "client")
        # self._key = Config().rcf.get("toollab", "authkey")
        # self._binds_id = Config().rcf.get("toollab", "binds_id")
        # self._interface_id = Config().rcf.get("toollab", "interface_id")
        # self._env_name = Config().rcf.get("toollab", "env_name")
        self._url = "http://apicenter.nsg.com/index/in"
        self._key = "458b97de4c0bb5bf416c8cea208309ed"
        self._binds_id = "5e8c0a091b1800007c006b1a"
        self._interface_id = 1348
        self._env_name = "offline"
