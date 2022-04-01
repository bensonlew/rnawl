# -*- coding: utf-8 -*-
# __author__ = 'wq'
# last modified @ 20200414
import json


def params_check(toollab_params):
    if not toollab_params["input_table"]:
        info = {"success": False, "info": '必须提供数据表'}
        return json.dumps(info)
    if not toollab_params["fq"]:
        info = {"success": False, "info": '必须设置频率'}
        return json.dumps(info)
    return None
