# -*- coding: utf-8 -*-
# __author__ = 'xuting'
import re
import json
import datetime
from collections import OrderedDict


def param_pack(param):
    if not isinstance(param, dict):
        raise Exception("传入的param不是一个字典")
    new_param = OrderedDict(sorted(param.items()))
    params = json.dumps(new_param)
    params = re.sub(':\s+', ':', params)
    params = re.sub(',\s+', ',', params)
    return params


def filter_json_sort(filter_detail):
    filters = json.loads(filter_detail)
    temp = []
    for i in filters:
        temp.append(OrderedDict(sorted(i.items(), key=lambda t: t[0])))
    return temp


def sub_group_detail_sort(detail):
    table_list = json.loads(detail)
    result_list = []
    for table_dict in table_list:
        if not isinstance(table_dict, dict):
            raise Exception("传入的table_dict不是一个字典")
        for keys in table_dict.keys():
            table_dict[keys] = sorted(table_dict[keys])
        sort_key = OrderedDict(sorted(table_dict.items(), key=lambda t: t[0]))
        result_list.append(sort_key)
    # result_list = json.dumps(result_list, sort_keys=True, separators=(',', ':'))
    return result_list


def group_detail_sort(detail):
    if isinstance(detail, dict):
        table_dict = detail
    else:
        table_dict = json.loads(detail)
    if not isinstance(table_dict, dict):
        raise Exception("传入的table_dict不是一个字典")
    for keys in table_dict.keys():
        table_dict[keys] = sorted(table_dict[keys])
    sort_key = OrderedDict(sorted(table_dict.items(), key=lambda t: t[0]))
    table_dict = sort_key
    # table_dict = json.dumps(table_dict, sort_keys=True, separators=(',', ':'))
    return table_dict


def get_lefse_catecory_name(detail):
    table_list = json.loads(detail)
    for table_dict in table_list:
        if not isinstance(table_dict, dict):
            raise Exception("传入的table_dict不是一个字典")
    if len(table_list) == 1:
        groupname = table_list[0].keys()
        groupname.sort()
        category = ','.join(groupname)
        second_category = ''
        return category, second_category
    else:
        groupname = table_list[0].keys()
        groupname.sort()
        category = ','.join(groupname)
        second_groupname = table_list[1].keys()
        second_groupname.sort()
        second_category = ','.join(second_groupname)
        return category, second_category


def GetUploadInfo(client, member_id, project_sn, task_id, name):
    if client == "client01":
        head = "sanger:"
        update_api = "meta.update_status"
    elif client == "client03":
        head = "tsanger:"
        update_api = "meta.tupdate_status"
    else:
        raise Exception("未知用户:{}".format(client))
    # strTime = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")  # modified by hesheng 20161223 名称与主表名称统一
    fullPath = "{}rerewrweset/files/{}/{}/{}/report_results/{}".format(head, member_id, project_sn, task_id, name)
    return (fullPath, update_api)


def GetUploadInfo_denovo(client, member_id, project_sn, task_id, name):
    if client == "client01":
        head = "sanger:"
        update_api = "denovo.update_status"
    elif client == "client03":
        head = "tsanger:"
        update_api = "denovo.tupdate_status"
    else:
        raise Exception("未知用户:{}".format(client))
    strTime = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    fullPath = "{}rerewrweset/files/{}/{}/{}/report_results/{}_{}".format(head, member_id, project_sn, task_id, name, strTime)
    return (fullPath, update_api)

def GetUploadInfo_pt(client, member_id, project_sn, task_id, name):
    if client == "client01":
        head = "sanger:"
        update_api = "paternity_test.update_status"
    elif client == "client03":
        head = "tsanger:"
        update_api = "paternity_test.tupdate_status"
    else:
        raise Exception("未知用户:{}".format(client))
    strTime = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    fullPath = "{}rerewrweset/files/{}/{}/{}/report_results/{}_{}".format(head, member_id, project_sn, task_id, name, strTime)
    return (fullPath, update_api)