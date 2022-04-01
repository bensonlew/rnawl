# -*- coding: utf-8 -*-

import json


def params_check(toollab_params):
    for arg in ["table"]:
        if arg not in toollab_params:
            info = {"success": False, "info": "缺失参数%s"%arg}
            return info

        method =["mean",'sum', 'median','none']
        if 'method' in toollab_params:
            if toollab_params['method'] not in method:
                info = {"success": False, "info": "method not in %s"%(' '.join(method))}
                return info
    return None
