# -*- coding: utf-8 -*-

import json


def params_check(toollab_params):
    for arg in ["table", "correct","sample1","sample2"]:
        if arg not in toollab_params:
            info = {"success": False, "info": "缺失参数%s"%arg}
            return info
    return None
