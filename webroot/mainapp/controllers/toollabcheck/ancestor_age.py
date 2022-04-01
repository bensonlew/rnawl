# -*- coding: utf-8 -*-

import json


def params_check(toollab_params):
    for arg in ["sample1","outgroup1","csample1","num_c_samples", "haplogroup"]:
        if arg not in toollab_params:
            info = {"success": False, "info": "缺失参数%s"%arg}
            return info