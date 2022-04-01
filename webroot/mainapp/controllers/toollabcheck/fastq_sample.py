# -*- coding: utf-8 -*-

import json


def params_check(toollab_params):
    if 'extract_type' not in toollab_params:
        return  {"success": False, "info": "缺失参数{}".format(a)}
    ex_t = toollab_params['extract_type'] + '_value'
    if ex_t not in toollab_params:
        return {"success": False, "info": "extract_type为{} 必须设置参数 {}".format(toollab_params['extract_type'], ex_t)}
    b1 = 'read1' in toollab_params
    b2 = 'read2' in toollab_params
    b = 'read' in toollab_params
    if not (b or b1 or b2):
        return {"success": False, "info": "缺失参数read1, read2 或read"}
    elif not (b1 and b2) and (b1 or b2):
        return {"success": False, "info": "缺失参数read1, read2必须成对配置"}
    return None
