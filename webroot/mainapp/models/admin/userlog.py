# -*- coding: utf-8 -*-
# __author__ = 'guoquan'  
from .base import BaseModel
import web


class UserlogModel(BaseModel):

    def __init__(self):
        super(UserlogModel, self).__init__()
        self._user = web.config.get('_session').user

    def add(self, info):
        return self.db.insert(self.table, user=self._user, discription=info)
