# -*- coding: utf-8 -*-

import json


def params_check(toollab_params):
    for arg in ["co_type", "co_value", "p_value", "top_abundance", "strategy"]:
        if arg not in toollab_params:
            info = {"success": False, "info": "缺失参数%s"%arg}
            return info
    if toollab_params['co_type'] not in ["pearson", "spearman", 'kendall']:
        return {'success': False, 'info': '不支持的相关系数类型{}, 目前只支持{}'.format(toollab_params['co_type'], ["pearson", "spearman", 'kendall'])}
    return None
