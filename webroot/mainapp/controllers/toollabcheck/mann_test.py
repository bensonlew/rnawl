# -*- coding: utf-8 -*-
# __author__ = 'wq'
# last modified @ 20200414
import json


def params_check(toollab_params):
    #if not toollab_params['mann_input']:
    #    info = {"success": False, "info": '必须提供数据表'}
    #    return json.dumps(info)
    #if not toollab_params['mann_group']:
    #    info = {"success": False, "info": '必须提供分组文件'}
    #    return json.dumps(info)
    if not toollab_params['sample']:
        info = {"success": False, "info": '必须设置样本名为行标签或列标签'}
        return json.dumps(info)
    if not toollab_params['compare_type']:
        info = {"success": False, "info": '必须设置秩和检验比较策略'}
        return json.dumps(info)
    else:
        if toollab_params['compare_type'] not in ["multi"]:
            if not toollab_params['group_name1']:
                info = {"success": False, "info": '必须设置样本组1名称'}
                return json.dumps(info)
            if not toollab_params['group_name2']:
                info = {"success": False, "info": '必须设置样本组2名称'}
                return json.dumps(info)
        else:
            if not toollab_params['group_name']:
                info = {"success": False, "info": '必须设置样本组名称'}
                return json.dumps(info)
    if toollab_params['correction'] not in ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"]:
        info = {"success": False, "info": '该多重检验校正的方法不被支持'}
        return json.dumps(info)
    # if toollab_params["coverage"]:
    #     if toollab_params["coverage"] not in [0.90, 0.95, 0.98, 0.99, 0.999]:
    #         info = {"success": False, "info": '秩和检验的置信区间的置信度不在范围值内'}
    #         return json.dumps(info)
    return None
