# -*- coding: utf-8 -*-
# __author__ = 'moli.zhou'
from .tupdate_status import TupdateStatus


class UpdateStatus(TupdateStatus):

    def __init__(self, data):
        super(UpdateStatus, self).__init__(data)
        self._client = "client01"
        # self._key = "1ZYw71APsQ"
        # self._url = "http://api.sanger.com/task/add_file"
        # self.database = self._mongo_client[self.config.MONGODB + '_paternity_test']
        self._binds_id = "5f4324a19b79000093008059"
        self._interface_id = 66
        self._env_name = "online"
        self._key = "e0bb04dd111155b2e6bc6db26d0e1fef"
        self._url = "http://apicenter.lab.majorbio.com/index/in"
