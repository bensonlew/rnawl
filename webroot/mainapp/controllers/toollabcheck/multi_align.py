# -*- coding: utf-8 -*-

import json


def params_check(toollab_params):
    for arg in ["fasta","align_method",'out_format']:
        if arg not in toollab_params:
            info = {"success": False, "info": "缺失参数%s"%arg}
            return info
    return None
