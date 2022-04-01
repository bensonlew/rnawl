# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# last modified @ 20200416
import json


def params_check(toollab_params):
    if not toollab_params["snp_table"]:
        info = {"success": False, "info": "必须输入snp_table文件"}
        return json.dumps(info)
    if not toollab_params["guide"]:
        info = {"success": False, "info": "必须输入是否显示显著性参考线"}
        return json.dumps(info)
    if not toollab_params["color1"]:
        info = {"success": False, "info": "必须输入颜色1"}
        return json.dumps(info)
    if not toollab_params["color2"]:
        info = {"success": False, "info": "必须输入颜色2"}
        return json.dumps(info)
    return None
