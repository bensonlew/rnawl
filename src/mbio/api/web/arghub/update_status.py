# -*- coding: utf-8 -*-
from ..meta.update_status import UpdateStatus
import os


class UpdateStatus(UpdateStatus):

    def __init__(self, data):
        super(UpdateStatus, self).__init__(data)
        self._project_type = 'arghub'


    def update_status(self):
        self.logger.info("###### update status####")
        self.logger.info(self.db)
        super(UpdateStatus, self).update_status()
        
