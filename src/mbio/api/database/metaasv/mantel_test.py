# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
from mainapp.libs.param_pack import group_detail_sort
import types

class MantelTest(Base):
    def __init__(self, bind_object):
        super(MantelTest, self).__init__(bind_object)
        self._project_type = 'metaasv'

    @report_check
    def add_mantel_table(self, level, otu_id, env_id, task_id=None, name=None, params=None, spname_spid=None):
        if level not in range(1, 10):
            self.bind_object.logger.error("level参数%s为不在允许范围内!" % level)
            self.bind_object.set_error("level参数不在范围内")
        origin_name = "mantel_origin"
        matrix_types = ["species_matrix", "env_matrix"]
        if "units" in params:
            matrix_types = ["species_matrix", "env_matrix", "partial_matrix"]
            origin_name = "partial_mantel_origin"

        if not isinstance(otu_id, ObjectId):
            otu_id = ObjectId(otu_id)
        collection = self.db["asv"]
        result = collection.find_one({"_id": otu_id})
        if task_id is None:
            task_id = result['task_id']
        if spname_spid:
            group_detail = {'All': [str(i) for i in spname_spid.values()]}
            params['group_detail'] = group_detail_sort(group_detail)
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": task_id,
            "env_id": env_id,
            "asv_id": otu_id,
            "name": self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else origin_name,
            "level_id": level,
            "status": "end",
            "desc": "",
            "matrix_types": matrix_types,
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["mantel"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_mantel_detail(self, file_path, mantel_id=None, main_colletion_name=None):
        """
        导入mantel的详情表
        :param file_path: 导表的数据
        :param mantel_id: 主表id
        :param main_colletion_name:
        :return:
        """
        if not isinstance(mantel_id, ObjectId):  # 检查传入的network_corfd_id是否符合ObjectId类型
            if isinstance(mantel_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                mantel_id = ObjectId(mantel_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('mantel_id必须为ObjectId对象或其对应的字符串！')
        main_colletion_time = ""
        if main_colletion_name:
            main_colletion_time = main_colletion_name[11:]
        data_list = []
        is_partial = 0
        with open(file_path, "r") as f:
            for line in f:
                if line.startswith("#"):continue
                elif line.startswith("DM1"):continue
                else:
                    line = line.strip().split("\t")
                    data = {
                        "mantel_id": mantel_id,
                        "dm1": "species_matrix"+main_colletion_time,
                        "dm2": "env_matrix"+main_colletion_time
                    }
                    if len(line) == 8:
                        data["entries_num"] = line[3]
                        data["permutation"] = line[6]
                        data["mantel_r"] = line[4]
                        data["p_value"] = line[5]
                        data["tail_type"] = line[7]
                        data["cdm"] = "partial_matrix"+main_colletion_time
                        is_partial = 1
                    else:
                        data["entries_num"] = line[2]
                        data["permutation"] = line[5]
                        data["mantel_r"] = line[3]
                        data["p_value"] = line[4]
                        data["tail_type"] = line[6]
                    data_list.append(data)
        try:
            collection = self.db["mantel_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入mantel检验结果数据出错:%s" % e)
            self.bind_object.set_error("导入mantel检验结果数据出错")
        else:
            self.bind_object.logger.info("导入mantel检验结果数据成功")
        try:
            main_collection = self.db["mantel"]
            if is_partial == 1:
                table_data = {"table_data": ["dm1", "dm2", "cdm", "mantel_r", "p_value", "permutation", "tail_type"]}
            else:
                table_data = {"table_data": ["dm1", "dm2", "mantel_r", "p_value", "permutation", "tail_type"]}
            table_data_json = json.dumps(table_data, sort_keys=True, separators=(',', ':'))
            settled_params = {"software" : "python-2.7"}
            settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
            main_collection.update({"_id": mantel_id}, {"$set":{"main_id": mantel_id,
                                                                "settled_params": settled_params_json,
                                                                "table_data": table_data_json}})
        except Exception, e:
            self.bind_object.logger.error("更新mantel检验结果数据出错:%s" % e)
            self.bind_object.set_error("更新mantel检验结果数据出错")
        else:
            self.bind_object.logger.info("更新mantel检验结果数据成功")


