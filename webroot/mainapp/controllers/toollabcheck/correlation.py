# -*- coding: utf-8 -*-
# __author__ = 'HD'
# last modified @ 20200407
import json


def params_check(toollab_params):
    if toollab_params["strategy"] not in ["sample", "group"]:
        info = {"success": False, "info": "strategy参数{}必须为sample or group".format(toollab_params["strategy"])}
        return json.dumps(info)
    if toollab_params["strategy"] == "group":
        if "group_method" not in toollab_params.keys():
            info = {"success": False, "info": "strategy 为 group的时候，必须要设置group_method参数！"}
            return json.dumps(info)
        if toollab_params['group_method'] not in ["mean", "sum", "median"]:
            info = {"success": False, "info": "group_method分组数据计算方法必须为mean or sum or median中任一个！"}
            return json.dumps(info)
    if toollab_params["distance_method"] not in ["euclidean", "maximum", "manhattan", "canberra", "binary",
                                                 "minkowski"]:
        info = {"success": False, "info": "不支持该距离算法"}
        return json.dumps(info)
    if toollab_params["method"] not in ["pearson", "spearman"]:
        info = {"success": False, "info": "不支持该相关系数方法"}
        return json.dumps(info)
    if toollab_params["cluster"] not in ["complete", "single", "average", "no"]:
        info = {"success": False, "info": "不支持该聚类方式"}
        return json.dumps(info)
    return None
