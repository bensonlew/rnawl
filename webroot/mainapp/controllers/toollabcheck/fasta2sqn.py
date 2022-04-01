# -*- coding: utf-8 -*-

import json


def params_check(toollab_params):
    for arg in ["temp_file", "seq_file"]:
        if arg not in toollab_params:
            info = {"success": False, "info": "缺失参数%s"%arg}
            return info
    a = int('structure_type' in toollab_params)
    b = int('structure_file' in toollab_params)
    if a + b == 1:
        return {'success': False, 'info': '必须同时设置structure_type和structure_file'}
    return None
