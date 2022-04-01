# -*- coding: utf-8 -*-
# __author__ = 'wq'
# last modified @ 20200413
import json


def params_check(toollab_params):
    if not toollab_params['tooltable']:
        info = {"success": False, "info": '必须提供数据表'}
        return json.dumps(info)
    return None
