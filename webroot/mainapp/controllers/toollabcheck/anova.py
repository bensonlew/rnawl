# -*- coding: utf-8 -*-

import json


def params_check(toollab_params):
    for arg in ["group_name","mul_test"]: #"table", "group_file",
        if arg not in toollab_params:
            info = {"success": False, "info": "缺失参数%s"%arg}
            return info

        # if '|' not in toollab_params['group_name'] and ',' not in toollab_params['group_name'] and '，' not in toollab_params['group_name']:
        #     info = {"success": False, "info": "group_name 值用 ,或|连接"}
        #     return info

        mul_test =["holm", "hochberg", "hommel", "bonferroni", "BH", "BY","FDR", "none"]
        if toollab_params['mul_test'] not in mul_test:
            info = {"success": False, "info": "mul_test not in %s"%(' '.join(mul_test))}
            return info

    return None
