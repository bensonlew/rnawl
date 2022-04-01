# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from ..meta.update_status import UpdateStatus


class DenovoUpdateStatus(UpdateStatus):

    def __init__(self, data):
        super(DenovoUpdateStatus, self).__init__(data)
        self.mongodb = self._mongo_client[self.config.MONGODB + '_rna']
