# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.api.database.base import Base, report_check
import re
import datetime
from random import choice
import json
import string
from bson.objectid import ObjectId
from types import StringTypes
from mainapp.libs.param_pack import group_detail_sort, param_pack
from mbio.packages.metaasv.common_function import find_group_name
from bson.son import SON


class PanCore(Base):
    def __init__(self, bind_object):
        super(PanCore, self).__init__(bind_object)
        self._project_type = 'metaasv'
        self.main_task_id = "_".join(self.bind_object.sheet.id.split('_')[0:2])
        self.unique_id = self.get_unique()

    @report_check
    def add_pan_core(self, params, level_id, asv_id, group_id=None, name=None, status=None, spname_spid=None):
        if asv_id and not isinstance(asv_id, ObjectId):
            if isinstance(asv_id, StringTypes):
                origin_asv_id = ObjectId(asv_id)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!")
        if group_id not in ["all", "All", "ALL"]:
            if not isinstance(group_id, ObjectId):
                if isinstance(group_id, StringTypes):
                    group_id = ObjectId(group_id)
                else:
                    self.bind_object.set_error("group_id必须为ObjectId对象或其对应的字符串!")
        collection = self.db["asv"]
        result = collection.find_one({"_id": origin_asv_id})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}在sg_otu表里找到相应的记录".format(str(asv_id)))
            self.bind_object.set_error("sg_otu表找不到相应记录")
        project_sn = result['project_sn']
        task_id = result['task_id']
        if spname_spid and params:
            if group_id not in [None, "All", "all", "ALL"]:
                ## 调用common模块，功能将导入的分组方案返回group_detail
                group_detail = find_group_name(task_id)
            else:
                group_detail = {'All': [str(i) for i in spname_spid.values()]}
            params['group_detail'] = group_detail_sort(group_detail)
        params = param_pack(params)
        if self.bind_object.sheet.main_table_name:
            name = self.bind_object.sheet.main_table_name
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "level_id": level_id,
            "asv_id": asv_id,
            "group_id": group_id,
            "status": "end",
            "desc": "PanCore_Origin",
            "unique_id": self.unique_id,
            "submit_location": "pancore",
            "name": name if name else "Pancore_Origin",
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["pan_core"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_pan_core_detail(self, file_path, pan_core_id, type):
        """
        导入详情表
        :param file_path: core/pan 的文件路径
        :param pan_core_id: 主表id string类型
        :param type: 导入表的类型
        :return:
        """
        if not isinstance(pan_core_id, ObjectId):
            if isinstance(pan_core_id, StringTypes):
                pan_core_id = ObjectId(pan_core_id)
            else:
                self.bind_object.set_error("pan_core_id必须为ObjectId对象或其对应的字符串!")
        with open(file_path, 'rb') as r:
            header = r.next().rstrip("\n")
            header = re.split('\t', header)
            header_list = header ##导入MongoDB的key
            pan_core_table = []
            pan_core_detail = []
            for line in r:
                line = line.rstrip('\n')
                line = re.split('\t', line)
                length = len(line)

                # 应前端要求，对category_name进行修改， 从group表的All改为all
                if line[0] in ["All", "all", "ALL"]:
                    temp_category_name = "all"
                else:
                    temp_category_name = line[0]
                insert_data = {
                    "pan_core_id": pan_core_id,
                    "specimen": temp_category_name,
                    "type" : type
                }
                for i in range(1, length):
                    insert_data[header_list[i]] = line[i]
                all_insert_data = SON(insert_data)
                pan_core_table.append(all_insert_data)
                for i in range(1, length):
                    data = {
                        "pan_core_id": pan_core_id,
                        "specimen": temp_category_name,
                        "type" : type,
                        "name" : ""
                    }
                    try:
                        data["x"] = int(header_list[i])
                    except:
                        data["x"] = 0
                    try:
                        data["y"] = float(line[i])
                        # if data["y"] != 0.0: ## 因为作图需要用到数据为0.0的数据，不然很难看--产品--张俊彪
                        all_data = SON(data)
                        pan_core_detail.append(all_data)
                    except:
                        if line[i] in ["NA", "na", "Na"]:
                            pass

            try:
                collection_table = self.db['pan_core_table']
                collection_table.insert_many(pan_core_table)
                self.bind_object.logger.info("导入pan_core_table表格成功！")
                collection_detail = self.db['pan_core_detail']
                collection_detail.insert_many(pan_core_detail)
                self.bind_object.logger.info("导入pan_core_detail表格成功！")
                if type in ["pan"]:
                    main_collection = self.db["pan_core"]
                    pan_core_type = ["pan", "core"]
                    settled_params = {"software" : "R-3.3.1 (vegan)"}
                    table_data = {
                        "table_data": ["specimen"] + header_list[1:],
                        "condition": {"type": pan_core_type}
                    }
                    table_data_json = json.dumps(table_data, sort_keys=True, separators=(',', ':'))
                    line_data = {
                        "line_data": {"name":"specimen",
                        "condition": {"type": pan_core_type}}
                    }
                    line_data_json = json.dumps(line_data, sort_keys=True, separators=(',', ':'))
                    main_collection.update({"_id":pan_core_id},{"$set":{"main_id": pan_core_id,
                                                                        "table_data": table_data_json,
                                                                        "line_data": line_data_json,
                                                                        "settled_params": json.dumps(settled_params)}})
                    self.bind_object.logger.info("更新pan_core表格成功！")

            except Exception as e:
                self.bind_object.logger.error("导入pan_core_detail表格{}信息出错:{}".format(file_path, e))
                self.bind_object.set_error("导入pan_cor_detail表出错")
            else:
                self.bind_object.logger.info("导入pan_core_detail表格{}成功".format(file_path))

    def get_unique(self):
        """
        get unique_id
        :return:
        """
        chars = string.ascii_letters + string.digits
        unique_id = "".join([choice(chars) for i in range(10)])
        collection = self.db["pan_core"]
        result = collection.find_one({"unique_id": unique_id, "task_id": self.main_task_id})
        if result:
            return self.get_unique()
        return unique_id
