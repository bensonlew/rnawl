# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
from .admin import Admin


class MainAction(Admin):

    def __init__(self):
        super(MainAction, self).__init__()
        self.nav_index = 0

    def GET(self):
        return self.render.admin.run(self)
