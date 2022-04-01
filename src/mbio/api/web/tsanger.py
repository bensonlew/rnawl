# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
import json
import urllib
from biocluster.wpm.log import Log
from biocluster.core.function import CJsonEncoder, filter_error_info
from .sanger import Sanger


class Tsanger(Sanger):

    def __init__(self, data):
        super(Tsanger, self).__init__(data)
        self._client = "client03"
        self._key = "458b97de4c0bb5bf416c8cea208309ed"
        self._binds_id = "5f4324f09b7900009300805c"
        self._interface_id = 57
        self._env_name = "offline"
        self._url = "http://apicenter.nsg.com/index/in"
        # self._url = "http://api.tsg.com/task/add_task_log"
        # self._post_data = "%s&%s" % (self.get_sig(), self.post_data)

    def update(self):
        self.send()
