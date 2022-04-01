# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.api.database.base import Base, report_check
import re
from bson.objectid import ObjectId
import datetime
import json
from bson.son import SON
from types import StringTypes
from mainapp.libs.param_pack import group_detail_sort, param_pack
from mbio.packages.metaasv.common_function import find_group_name


class Estimator(Base):
    def __init__(self, bind_object):
        super(Estimator, self).__init__(bind_object)
        self._project_type = 'metaasv'
        self.main_task_id = "_".join(self.bind_object.sheet.id.split('_')[0:2])

    @report_check
    def add_est_table(self, file_path, level, major=False, otu_id=None, est_id=None, task_id=None, name=None,
                      params=None, spname_spid=None, indices=None,group_id=None):
        if level not in range(1, 10):
            self.bind_object.logger.error("level参数%s为不在允许范围内!" % level)
            self.bind_object.set_error("level不在允许范围内")
        data_list = []
        if not isinstance(otu_id, ObjectId):
            otu_id = ObjectId(otu_id)
        task_id = self.main_task_id
        if major:
            params['asv_id'] = str(otu_id)
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
                "name": name if name else "Estimators_Origin",
                "level_id": level,
                "status": "end",
                "desc": "",
                "params": params,
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
            collection = self.db["alpha_diversity"]
            est_id = collection.insert_one(insert_data).inserted_id
        else:
            if est_id is None:
                self.bind_object.set_error("major为False时需提供est_id!")
            if not isinstance(est_id, ObjectId):
                if isinstance(est_id, StringTypes):
                    est_id = ObjectId(est_id)
                else:
                    self.bind_object.set_error("est_id必须为ObjectId对象或其对应的字符串!")
        with open(file_path, 'r') as f:
            l = f.readline().strip('\n')
            if not re.match(r"^sample", l):
                self.bind_object.logger.error("文件%s格式不正确，请选择正确的estimator表格文件" % file_path)
                self.bind_object.set_error("estimator表格格式不正确")
            est_list = l.split("\t")
            est_list.pop(0)
            while True:
                line = f.readline().strip('\n')
                if not line:
                    break
                line_data = line.split("\t")
                sample_name = line_data.pop(0)
                data = [("alpha_diversity_id", est_id), ("specimen", sample_name)]
                i = 0
                for est in est_list:
                    data.append((est, float(line_data[i])))
                    i += 1
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["alpha_diversity_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["alpha_diversity"]
            settled_params = {'software': "mothur-1.30"}
            settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
            index_list = indices.split(",")
            table_data = {"table_data": ["specimen"] + index_list}
            table_data_json = json.dumps(table_data, sort_keys=True, separators=(',', ':'))
            column_data = {
                "column_data": {"name":"specimen",
                        "data":"mean",
                        "condition": {"type":index_list}
                }}
            column_data_json = json.dumps(column_data, sort_keys=True, separators=(',', ':'))
            ishape_data = {
                "ishape_data": {"name":"specimen",
                        "data":["mean", "lower", "upper"],
                        "condition": {"type":index_list}
                }}
            ishape_data_json = json.dumps(ishape_data, sort_keys=True, separators=(',', ':'))
            main_collection.update({"_id":est_id},{"$set":{"main_id": est_id,
                                                        "table_data": table_data_json,
                                                        "column_data": column_data_json,
                                                        "ishape_data": ishape_data_json,
                                                        "settled_params": settled_params_json}})
        except Exception, e:
            self.bind_object.logger.error("导入estimator%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入estimator%s信息成功!" % file_path)
        return est_id

    def add_est_bar(self, file_path, est_id, indices=None):
        data_list = []
        if not isinstance(est_id, ObjectId):
            if isinstance(est_id, StringTypes):
                est_id = ObjectId(est_id)
            else:
                self.bind_object.set_error("est_id必须为ObjectId对象或其对应的字符串!")
        with open(file_path, 'r') as f:
            l = f.readline().strip('\n')
            if not re.match(r"^sample", l):
                self.bind_object.logger.error("文件%s格式不正确，请选择正确的estimator表格文件" % file_path)
                self.bind_object.set_error("estimator表格格式不正确")
            est_list = l.split("\t")[1:]
            for line in f:
                line_data = line.strip().split("\t")
                sample_name = line_data[0]
                data = [("alpha_diversity_id", est_id), ("specimen", sample_name)]
                indices_list = indices.split(",")
                for indic in indices_list:
                    if indic in est_list:
                        indic_index = est_list.index(indic)
                        data.append(("mean", float(line_data[indic_index+1])))
                    indic_lci = indic + "_lci"
                    if indic_lci in est_list:
                        indic_index = est_list.index(indic_lci)
                        data.append(("lower", float(line_data[indic_index+1])))
                    indic_hci = indic + "_hci"
                    if indic_hci in est_list:
                        indic_index = est_list.index(indic_hci)
                        data.append(("upper", float(line_data[indic_index+1])))
                    data.append(("type", indic))
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            collection = self.db["alpha_diversity_bar"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入estimator%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入estimator%s信息成功!" % file_path)
        return est_id