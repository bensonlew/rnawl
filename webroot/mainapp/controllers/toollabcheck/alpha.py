# -*- coding: utf-8 -*-
# __author__ = 'zzg'
# last modified @ 20210316
import json


def params_check(toollab_params):
    if 'data_table' not in toollab_params:
        info = {"success": False, "info": "参数{}不存在".format('data_table')}
        return json.dumps(info)
    return None
