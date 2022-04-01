# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
import web
import os
from mainapp.models.admin.userlog import UserlogModel


class LogoutAction(object):

    def GET(self):
        session = web.config.get('_session')
        if session.user:
            model = UserlogModel()
            model.add("退出登陆！")
        session.user = None
        session.is_login = False
        session.is_admin = False
        path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../views/'))
        render = web.template.render(path)
        return render.admin.logout()
