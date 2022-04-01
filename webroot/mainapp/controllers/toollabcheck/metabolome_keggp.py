# -*- coding: utf-8 -*-

import json


def params_check(toollab_params):
    for arg in ["metab_name"]:
        if arg not in toollab_params:
            info = {"success": False, "info": "缺失参数%s"%arg}
            return info
    if toollab_params['metab_name'].strip() == '':
        info = {"success": False, "info": "输入数据不能为空"}
        return info
    return None
