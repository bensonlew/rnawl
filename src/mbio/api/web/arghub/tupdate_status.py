# -*- coding: utf-8 -*-
from .update_status import UpdateStatus


class TupdateStatus(UpdateStatus):

    def __init__(self, data):
        super(TupdateStatus, self).__init__(data)
        self._client = "client03"
        self._project_type = 'arghub'
        self._env_name = "offline"
        self._binds_id = '60890888f766300c19191ac1'
        self._key = "458b97de4c0bb5bf416c8cea208309ed"
        self._url = "http://apicenter.nsg.com/index/in"
