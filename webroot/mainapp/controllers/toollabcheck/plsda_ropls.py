# -*- coding: utf-8 -*-

import json


def params_check(toollab_params):
    for arg in ["table","group_table","method"]:
        if arg not in toollab_params:
            info = {"success": False, "info": "缺失参数%s"%arg}
            return info
    if toollab_params['method'] not in ['plsda', 'oplsda']:
        info = {"success": False, "info": "method not in [plsda, oplsda]"}
        return info
    return None
