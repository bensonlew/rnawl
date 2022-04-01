# -*- coding: utf-8 -*-
# __author__ = 'guoquan'  
from .admin import Admin
import web


class IndexAction(Admin):

    def __init__(self):
        super(IndexAction, self).__init__()

    def GET(self):
        raise web.seeother('main')
