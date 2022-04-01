# -*- coding: utf-8 -*-
# __author__ = 'xuting'
# last modified guhaidong 20171116

from biocluster.api.database.base import Base, report_check
import re
import datetime
from random import choice
import json
import string
from bson.objectid import ObjectId
from types import StringTypes
# from biocluster.config import Config
from mainapp.libs.param_pack import group_detail_sort, param_pack


class PanCore(Base):
    def __init__(self, bind_object):
        super(PanCore, self).__init__(bind_object)
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB
        self.main_task_id = "_".join(self.bind_object.sheet.id.split('_')[0:2])  # add main_task_id by guhaidong 20171116
        self.unique_id = self.get_unique()

    @report_check
    def create_pan_core_table(self, pan_core_type, params, group_id, level_id, from_otu_table=0, name=None, status=None, spname_spid=None):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!", code="51004801")
        if group_id not in ["all", "All", "ALL"]:
            if not isinstance(group_id, ObjectId):
                if isinstance(group_id, StringTypes):
                    group_id = ObjectId(group_id)
                else:
                    self.bind_object.set_error("group_id必须为ObjectId对象或其对应的字符串!", code="51004802")
        if not status:
            status = "end"
        if spname_spid:
            my_params = json.loads(params)
            group_detail = {'All': [str(i) for i in spname_spid.values()]}
            my_params['group_detail'] = group_detail_sort(group_detail)
            my_params['otu_id'] = str(from_otu_table)
            params = param_pack(my_params)
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": from_otu_table})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}在sg_otu表里找到相应的记录".format(str(from_otu_table)))
            self.bind_object.set_error("sg_otu表找不到相应记录", code="51004803")
        project_sn = result['project_sn']
        task_id = result['task_id']
        if self.bind_object.sheet.main_table_name:
            name = self.bind_object.sheet.main_table_name
            name = name.split('_')
            if pan_core_type == 1:
                name.pop(1)
            else:
                name.pop(0)
            name = '_'.join(name)
        if pan_core_type == 1:
            desc = "正在计算pan otu表格"
        else:
            desc = "正在计算core otu表格"
        insert_data = {
            "type": pan_core_type,
            "project_sn": project_sn,
            "task_id": task_id,
            "level_id": level_id,
            "otu_id": from_otu_table,
            "group_id": group_id,
            "status": status,
            "desc": desc,
            "unique_id": self.unique_id,
            "submit_location": "otu_pan_core",
            "name": name if name else "pan_core表格",
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_otu_pan_core"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_pan_core_detail(self, file_path, pan_core_id):
        if not isinstance(pan_core_id, ObjectId):
            if isinstance(pan_core_id, StringTypes):
                pan_core_id = ObjectId(pan_core_id)
            else:
                self.bind_object.set_error("pan_core_id必须为ObjectId对象或其对应的字符串!", code="51004804")
        with open(file_path, 'rb') as r:
            header = r.next().rstrip("\n")
            header = re.split('\t', header)
            pan_core_detail = list()
            for line in r:
                line = line.rstrip('\n')
                line = re.split('\t', line)
                value = list()
                length = len(line)
                for i in range(1, length):
                    my_dic = dict()
                    my_dic["title"] = header[i]
                    my_dic["value"] = line[i] if line[i] != "NA" else 0
                    value.append(my_dic)

                # 应前端要求，对category_name进行修改， 从group表的All改为all
                if line[0] in ["All", "all", "ALL"]:
                    temp_category_name = "all"
                else:
                    temp_category_name = line[0]
                insert_data = {
                    "pan_core_id": pan_core_id,
                    "category_name": temp_category_name,
                    "value": value
                }
                pan_core_detail.append(insert_data)
            try:
                collection = self.db['sg_otu_pan_core_detail']
                collection.insert_many(pan_core_detail)

                main_collection = self.db["sg_otu_pan_core"]
                #main_collection.update({"_id":pan_core_id},{"$set":{"main_id": pan_core_id}})

            except Exception as e:
                self.bind_object.logger.error("导入pan_core_detail表格{}信息出错:{}".format(file_path, e))
                self.bind_object.set_error("导入pan_cor_detail表出错", code="51004805")
            else:
                self.bind_object.logger.info("导入pan_core_detail表格{}成功".format(file_path))

    def get_unique(self):
        chars = string.ascii_letters + string.digits
        unique_id = "".join([choice(chars) for i in range(10)])
        collection = self.db["sg_otu_pan_core"]
        result = collection.find_one({"unique_id": unique_id, "task_id": self.main_task_id})
        if result:
            return self.get_unique()
        return unique_id
