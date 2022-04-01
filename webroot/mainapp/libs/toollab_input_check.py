# -*- coding: utf-8 -*-
# __author__ = 'HD'
# last modified @ 20200407
import functools
import json
import web


def check_format(f):
    @functools.wraps(f)
    def wrapper(obj):
        data = web.input()
        params = ["params", 'basis']
        for param in params:
            if not hasattr(data, param):
                info = {"success": False, "info": "缺少%s参数!" % param}
                return json.dumps(info)
        basis = json.loads(data.basis)
        anaysis_name = ["type", "task_type", "name", "main_table_name", "project_sn", "task_id"]
        for name in anaysis_name:
            if name not in basis.keys():
                info = {"success": False, "info": "basis参数中缺少%s子参数!" % name}
                return json.dumps(info)
        if basis["type"] not in ['tool', 'module', 'workflow']:
            info = {"success": False, "info": "type--{}，参数必须为tool，module， workflow".format(basis["type"])}
            return json.dumps(info)
        if basis["task_type"] not in ["instant", "submit"]:
            info = {"success": False, "info": "参数 %s 不合法, 必须为instant or submit！" % basis["task_type"]}
            return json.dumps(info)
        if not any(json.loads(data.params)):
            info = {"success": False, "info": "params对象中没有值，查看小工具后台参数配置是否正确！"}
            return json.dumps(info)
        return f(obj)
    return wrapper
