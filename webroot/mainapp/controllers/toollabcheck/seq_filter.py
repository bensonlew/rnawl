# -*- coding: utf-8 -*-
# __author__ = 'fwy'
import json
import re


def params_check(toollab_params):
    max_length = float(toollab_params["max_len"]) if "max_len" in toollab_params else 100000000
    min_length = float(toollab_params["min_len"]) if "min_len" in toollab_params else 0
    if max_length <= min_length:
        info = {'success': False,
                'info': "The maximum must be greater than the minimum!"}
        return json.dumps(info)
    return None
