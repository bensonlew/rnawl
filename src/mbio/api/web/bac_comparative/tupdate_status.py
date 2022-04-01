# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
from biocluster.config import Config
from ..meta.update_status import UpdateStatus


class TupdateStatus(UpdateStatus):

    def __init__(self, data):
        super(TupdateStatus, self).__init__(data)
        self._client = "client03"
        self._key = "hM4uZcGs9d"
        self._url = "http://api.tsg.com/task/add_file"
        self._project_type = 'bac_comparative'