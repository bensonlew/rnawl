# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
from ..meta.update_status import UpdateStatus


class UpdateStatus(UpdateStatus):

    def __init__(self, data):
        super(UpdateStatus, self).__init__(data)
        self._project_type = 'metagbin'
