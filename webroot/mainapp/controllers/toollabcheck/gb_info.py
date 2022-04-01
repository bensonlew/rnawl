# -*- coding: utf-8 -*-

import json


def params_check(toollab_params):
    for arg in ['gb_file', 'extract']:
        if arg not in toollab_params:
            info = {"success": False, "info": "缺失参数%s"%arg}
            return info
    a = set(toollab_params['extract'].split(',')) - set(['cds', 'prot', 'gff', 'genome'])
    if a:
        return {'success': False, 'info': '存在不支持的提取类型{}, 支持的提取类型为{}'.format(a, ['cds', 'prot', 'gff', 'genome'])}
    return None
