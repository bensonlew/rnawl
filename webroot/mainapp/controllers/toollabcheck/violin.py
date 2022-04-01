# -*- coding: utf-8 -*-
# __author__ = 'wq'
# last modified @ 20200514
import json


def params_check(toollab_params):
    if not toollab_params['tooltable']:
        info = {"success": False, "info": '必须提供数据表'}
        return json.dumps(info)
    # if not toollab_params['grouptable_checked']:
    #     info = {"success": False, "info": '必须选择是否上传分组文件'}
    #     return json.dumps(info)
    # if toollab_params['grouptable_checked']:
    #     if toollab_params['grouptable_checked'] in ["true"]:
    #         if not toollab_params['grouptable']:
    #             info = {"success": False, "info": '请上传分组文件'}
    #             return json.dumps(info)
    return None
