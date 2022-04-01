# -*- coding: utf-8 -*-
# __author__ = 'guoquan'  

import web
import os
from mainapp.models.admin.base import BaseListModel, BaseModel
from mainapp.libs.pagination import Paginatiton
import re
import importlib
import json
from biocluster.config import Config


path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../views/'))
render = web.template.render(path)


def nologin():
    _is_ajax = False
    _request_type = web.ctx.env.get('X-Requested-With')
    if _request_type == "XMLHttpRequest":
        _is_ajax = True
    if _is_ajax:
        return web.Unauthorized(json.dumps({"success": False, "info": "没有登陆或登陆超时！"}))
    else:
        return web.Unauthorized(render.admin.nologin())


web.unauthorized = nologin


def noright():
    _is_ajax = False
    _request_type = web.ctx.env.get('X-Requested-With')
    if _request_type == "XMLHttpRequest":
        _is_ajax = True
    if _is_ajax:
        return web.Unauthorized(json.dumps({"success": False, "info": "没有权限访问！"}))
    else:
        return web.Forbidden(render.admin.noright())


web.forbidden = noright


class Admin(object):
    def __init__(self):
        self.title = "后台管理"
        self.is_admin = 1 if self.session.is_admin else 0
        self.render = web.template.render(path, base='admin/layout')
        self.check_login()
        self.username = self.session.user
        self.nav_index = 0
        self._table = None
        self._pager = None
        self._list_model = None
        self._params = None
        self._error_render = render
        self.config = Config()
        self.info = ""
        self.url = ""

    @property
    def session(self):
        return web.config.get('_session')

    def log(self, info):
        model = self.get_model("userlog")
        model.add(info)

    def _get_table_name(self):
        class_name = self.__class__.__name__
        return re.sub(r"Action$", "", class_name).lower()

    @property
    def table(self):
        if not self._table:
            self._table = self._get_table_name()
        return self._table

    def check_admin(self):
        if not self.session.is_admin:
            self._set_noright()

    def check_login(self):
        if not self.session.is_login:
            self._set_unlogin()

    def _set_unlogin(self):
        if hasattr(self, "GET"):
            setattr(self, "GET", self._raise_nologin)
        if hasattr(self, "POST"):
            setattr(self, "POST", self._raise_nologin)

    def _raise_nologin(self):
        raise web.unauthorized()

    def _set_noright(self):
        if hasattr(self, "GET"):
            setattr(self, "GET", self._raise_noright)
        if hasattr(self, "POST"):
            setattr(self, "POST", self._raise_noright)

    def _raise_noright(self):
        raise web.forbidden()

    def _set_failed(self):
        if hasattr(self, "GET"):
            setattr(self, "GET", self._raise_failed)
        if hasattr(self, "POST"):
            setattr(self, "POST", self._raise_failed)

    def _raise_failed(self):
        return self.failed(self.info, self.url)

    @property
    def list_model(self):
        if not self._list_model:
            self._list_model = BaseListModel(self.table)
        return self._list_model

    @staticmethod
    def get_model(model_name):
        try:
            imp = importlib.import_module("mainapp.models.admin.%s" % model_name)
        except ImportError:
            return BaseModel(model_name)
        else:
            return getattr(imp, "%sModel" % model_name.capitalize())()

    @property
    def pager(self):
        if not self._pager:
            self._pager = Paginatiton(self.list_model)
        return self._pager

    @property
    def current_page(self):
        data = web.input()
        page = data.page if hasattr(data, "page") else 1
        try:
            page = int(page)
        except Exception:
            page = 1
        return page

    @property
    def params(self):
        if not self._params:
            self._params = Params(web.input())
        return self._params

    def check_required(self, required_list, relation="and"):
        if relation == "or":
            check_pass = False
            for p in required_list:
                if getattr(self.params, p) is not None and getattr(self.params, p) != "":
                    check_pass = True
            if not check_pass:
                self.info = "缺少参数！"
                return self._set_failed()
        else:
            for p in required_list:
                if getattr(self.params, p) is None or getattr(self.params, p) == "":
                    self.info = "缺少参数！"
                    return self._set_failed()

    def failed(self, info, url=None):
        _is_ajax = False
        _request_type = web.ctx.env.get('X-Requested-With')
        if _request_type == "XMLHttpRequest":
            _is_ajax = True
        if _is_ajax:
            return json.dumps({"success": False, "info": info})
        else:
            return self._error_render.admin.failed(info, url)


class Params(object):
    def __init__(self, data):
        self._data = data

    def __getattr__(self, name):
        if hasattr(self._data, name):
            return getattr(self._data, name)
        else:
            return None
