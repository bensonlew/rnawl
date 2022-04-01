# -*- coding: utf-8 -*-
# __author__ = 'wq'
# last modified @ 20200506
import json


def params_check(toollab_params):
    if not toollab_params['tooltable']:
        info = {"success": False, "info": '必须提供数据表'}
        return json.dumps(info)
    if not toollab_params['x_axis']:
        info = {"success": False, "info": '必须设置x轴数据'}
        return json.dumps(info)
    if not toollab_params['y_axis']:
        info = {"success": False, "info": '必须设置y轴数据'}
        return json.dumps(info)
    if not toollab_params['size']:
        info = {"success": False, "info": '必须设置气泡大小'}
        return json.dumps(info)
    return None
