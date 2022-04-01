# -*- coding: utf-8 -*-
# __author__ = 'zzg'
# last modified @ 20210316
import json


def params_check(toollab_params):
    if 'data_table' not in toollab_params:
        info = {"success": False, "info": "参数{}不存在".format('data_table')}
        return json.dumps(info)
    if 'in_fasta' not in toollab_params:
        info = {"success": False, "info": "参数{}不存在".format('data_table')}
        return json.dumps(info)
    if toollab_params["style"] not in ["16S", "18S","ITS"]:
        info = {"success": False, "info": "strategy参数{}必须为16S,18S或者ITS".format(toollab_params["style"])}
        return json.dumps(info)
    return None
