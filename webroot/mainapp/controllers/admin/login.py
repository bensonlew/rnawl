# -*- coding: utf-8 -*-
# __author__ = 'guoquan'  
import web
import os
from mainapp.models.user import User
from mainapp.models.admin.userlog import UserlogModel
import hashlib

session = web.config.get('_session')
path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../views/'))
render = web.template.render(path)


class LoginAction(object):

    def __init__(self):
        super(LoginAction, self).__init__()

    def GET(self):
        return render.admin.login()

    def POST(self):
        data = web.input()
        if not (hasattr(data, "username") and hasattr(data, "password")):
            raise web.unauthorized()
        user = User(data.username)
        m = hashlib.md5()
        m.update(data.password)
        if user.check_pass(m.hexdigest()):
            session.user = data.username
            session.is_login = True
            model = UserlogModel()
            model.add("登陆成功！")
            user.login()
            if user.is_admin():
                session.is_admin = True
            web.redirect("main")
        else:
            return render.admin.loginfailed()

