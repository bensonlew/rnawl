# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
from biocluster.api.database.base import Base, report_check
import os
import datetime
from bson.objectid import ObjectId
from types import StringTypes
import json
from mainapp.libs.param_pack import group_detail_sort,param_pack
from mbio.packages.metaasv.common_function import find_group_name
import copy
import numpy as np


class Rarefaction(Base):
    def __init__(self, bind_object):
        super(Rarefaction, self).__init__(bind_object)
        self._project_type = 'metaasv'
        self.category_x = []

    @report_check
    def add_rarefaction_detail(self, rare_id, file_path, indices,group=None):
        """
        根据需求，增加了group组内计算的方式
        :param rare_id:
        :param file_path:
        :param indices:
        :return:
        """
        if not isinstance(rare_id, ObjectId):
            if isinstance(rare_id, StringTypes):
                rare_id = ObjectId(rare_id)
            else:
                self.bind_object.set_error("rarefaction_id必须为ObjectId对象或其对应的字符串!")
        collection_first = self.db['rarefaction']
        result = collection_first.find_one({"_id": rare_id})
        if group:
            group_table = group
        else:
            group_table = self.bind_object.option("group").prop['path']

        # data_list = []
        sample_group_dict = {}##获取样本与group对应关系
        group_list = []##获取所有分组情况
        with open(group_table, 'r') as group_f:
            lins = group_f.readlines()
            for lin in lins[1:]:
                lin = lin.strip().split("\t")
                sample_na = lin[0]
                group_name = lin[1]
                sample_group_dict[sample_na] = group_name
                if group_name not in group_list:
                    group_list.append(group_name)

        task_id = result['task_id']
        rare_paths = os.listdir(file_path)
        rare_detail = []
        category_x = []
        sample_list = []
        all_group_data = []
        for rare_path in rare_paths:
            files = os.listdir(os.path.join(file_path, rare_path))
            fs_path = []
            all_x = []
            for f in files:
                fs_path.append(os.path.join(file_path, rare_path, f))
            all_group_sample = {}
            for fs in fs_path:
                sample_name = fs.split('.')[1]
                sample_name = sample_name.strip()
                # self.bind_object.logger.info("sample_name: {}".format(sample_name))
                if sample_name  not in sample_list:
                    sample_list.append(sample_name)
                if sample_name in sample_group_dict.keys():
                    group_name = sample_group_dict[sample_name]
                else:
                    self.bind_object.set_error("样本名称没有对应的分组名称，请检查group表！")
                x_y_dict = {}
                sample_x = {}
                # self.bind_object.logger.info("group_name: {}".format(group_name))
                with open(fs) as f:
                    lines = f.readlines()
                    for line in lines[1:]:
                        line_data = line.strip().split("\t")
                        if line_data[0] not in category_x:
                            category_x.append(line_data[0])
                        if line_data[0] not in all_x:
                            all_x.append(line_data[0])
                        x_y_dict[line_data[0]] = line_data[1]
                        insert_data = {
                            "rarefaction_id": rare_id,
                            "type": rare_path,
                            "specimen": sample_name,
                            "x": int(line_data[0]),
                            "y": float(line_data[1]),
                            "group": ""
                        }
                        rare_detail.append(insert_data)
                if group_name in all_group_sample.keys():
                    all_sample_list = all_group_sample[group_name]
                else:
                    all_sample_list = []
                if sample_x not in all_sample_list:
                    sample_x[sample_name] = x_y_dict
                    all_sample_list.append(sample_x)
                all_group_sample[group_name] = all_sample_list
            ##下面开始导组内合并的数据
            self.bind_object.logger.info("group_list:{}++++++++++++++".format(group_list))
            # self.bind_object.logger.info("sample_group_dict:{}+++++++++++++".format(sample_group_dict))
            # self.bind_object.logger.info("all_x:{}+++++++++++++++++++++".format(all_x))
            # self.bind_object.logger.info("all_group_sample:{}--------------------++++++++++++++".format(all_group_sample))
            # self.bind_object.logger.info("rare_path: {}-----------+++++++++++++".format(rare_path))
            data_list = self.calculate_group(group_list, all_x, all_group_sample, rare_path,rare_id)
            all_group_data += data_list
        try:
            collection = self.db['rarefaction_line']
            collection.insert_many(rare_detail)
            collection.insert_many(all_group_data)
            indice_list = indices.strip().split(",")
            settled_params = {"software" : "mothur-1.30"}
            settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
            line_data = {
                "line_data": {"name":"specimen",
                        "condition": {"type":indice_list}
                }
            }
            line_data_json = json.dumps(line_data, sort_keys=True, separators=(',', ':'))
            collection_first.update({"_id": ObjectId(rare_id), "task_id": task_id},
                                    {"$set": {"category_x": category_x,
                                              "main_id": ObjectId(rare_id),
                                              "line_data": line_data_json,
                                              "settled_params": settled_params_json,
                                              "specimen_list": sample_list,
                                              "status": "end"}})
        except Exception as e:
            self.bind_object.logger.error("导入rarefaction_line表格{}信息出错:{}".format(file_path, e))
            self.bind_object.set_error("导入rarefaction_line表格出错")
        else:
            self.bind_object.logger.info("导入rarefaction_line表格{}成功".format(file_path))


    def calculate_group(self, group_list, all_x, group_dict, indices, rare_id):
        """
        计算组内的合并方式，mean和median
        :param group_list: 组的list
        :param rare_id: 主表的_id
        :param all_x: 横坐标的list
        :param group_dict: 组、样本、横坐标对应关系
        :param indices: 指数
        :return:
        """
        data_list = []
        group_list = list(set(group_list))
        for group in group_list:
            x_y_dict = {}
            for x in all_x:
                sample_x_list = []
                group_sample_list = group_dict[group]
                for sample_dict in group_sample_list:
                    sample_name = sample_dict.keys()[0]
                    value_dict  = sample_dict[sample_name]
                    if x in value_dict.keys():
                        value = float(value_dict[x])
                        sample_x_list.append(value)
                x_y_dict[x] = sample_x_list
                    # else:
                    #     # self.bind_object.logger.info("x: {}-----------+++++++++++++".format(x))
                    #     all_group_list.append(0.0)
                # self.bind_object.logger.info("x: {}-----------+++++++++++++".format(all_group_list))
            self.bind_object.logger.info("all_group_list: {}-----------+++++++++++++".format(x_y_dict))
            if len(x_y_dict) != 0:
                for value_x in x_y_dict:
                    value_y = x_y_dict[value_x]
                    value_y = list(value_y)
                    if len(value_y) != 0:
                        average_value = np.mean(value_y)
                        median_value = np.median(value_y)
                        max_value = max(value_y)
                        min_value = min(value_y)
                        std_value = np.std(value_y)
                        insert_data = {
                            "x": int(value_x),
                            "mean": float(average_value),
                            "lower": float(average_value) - float(std_value),
                            "upper": float(average_value) + float(std_value),
                            "type": indices,
                            "group": "mean",
                            "rarefaction_id": rare_id,
                            "specimen": group
                        }
                        data_list.append(insert_data)
                        if (median_value == 0.0) and (int(value_x) == 1):
                            insert_data2 = {
                                "x": int(value_x),
                                "median": float(median_value),
                                "min": float(min_value),
                                "max": float(max_value),
                                "type": indices,
                                "group": "median",
                                "rarefaction_id": rare_id,
                                "specimen": group
                            }
                            data_list.append(insert_data2)
                        elif median_value != 0.0:
                            insert_data2 = {
                                    "x": int(value_x),
                                    "median": float(median_value),
                                    "min": float(min_value),
                                    "max": float(max_value),
                                    "type": indices,
                                    "group": "median",
                                    "rarefaction_id": rare_id,
                                    "specimen": group
                                }
                            data_list.append(insert_data2)
        return data_list

    @report_check
    def add_rare_table(self, file_path, level, otu_id=None, task_id=None, name=None, params=None, spname_spid=None,group_id=None):
        if level not in range(1, 10):
            self.bind_object.logger.error("level参数%s为不在允许范围内!" % level)
            self.bind_object.set_error("level参数不在允许的范围内")
        if task_id is None:
            task_id = self.bind_object.sheet.id
        if otu_id:
            if not isinstance(otu_id, ObjectId):
                otu_id = ObjectId(otu_id)
            params['asv_id'] = str(otu_id)  # otu_id在再metabase中不可用
        if spname_spid and params:
            if group_id not in [None, "All", "all", "ALL"]:
                ## 调用common模块，功能将导入的分组方案返回group_detail
                group_detail = find_group_name(task_id)
            else:
                group_detail = {'All': [str(i) for i in spname_spid.values()]}
            params['group_detail'] = group_detail_sort(group_detail)
        params = param_pack(params)
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": task_id,
            "asv_id": otu_id,
            "name": name if name else "Rarefaction_Origin",
            "level_id": level,
            "status": "start",
            "desc": "",
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["rarefaction"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        # self.add_rarefaction_detail(inserted_id, file_path)
        return inserted_id
