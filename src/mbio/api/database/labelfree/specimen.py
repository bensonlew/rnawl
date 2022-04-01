# -*- coding: utf-8 -*-
# __author__ = 'litangjian'
# modified by shicaiping 20180322

import json
import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.labelfree.api_base import ApiBase
import re
from collections import OrderedDict
import os
from bson.son import SON

# 数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到
class Specimen(ApiBase):
    def __init__(self, bind_object):
        super(Specimen, self).__init__(bind_object)
        self._project_type = 'labelfree'

    #@report_check
    def add_specimen(self, group_txt=None, params=None,
                project_sn='labelfree', task_id='labelfree', type='tmt'):
        if params is None:
            params = {"software": "peaks 8.5"}
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        data_list = []
        with open(group_txt, "r") as f:
            first_line = f.readline()
            if not re.match(r"#sample", first_line):
                raise Exception("%s文件类型不正确" % group_txt)
            for line in f:
                line = line.split()
                insert_data = {
                    'project_sn': project_sn,
                    'task_id': task_id,
                    'params': params if params else "",
                    'status': 'start',
                    'desc': '样品信息分组表',
                    'created_ts': datetime.datetime.now().strftime('%Y-%m-%d ''%H:%M:%S'),
                    "old_name": line[0],
                    "new_name": line[0],
                    "desc" : "",
                    "group" : line[1],
                    "type": type,
                }
                data_list.append(insert_data)
        try:
            collection = self.db["sg_specimen"]
            result = collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.set_error("导入样品信息数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入样品信息数据成功")
            self.sample_ids = result.inserted_ids
        sample_ids = [str(i) for i in result.inserted_ids]
        return sorted(sample_ids)

    def add_specimen_group(self, file):
        category_names = list()
        specimen_names = list()
        group_dict = OrderedDict()
        with open(file, "r") as f:
            f.readline()
            for line in f:
                tmp = line.strip().split("\t")
                group_dict.setdefault(tmp[1], list())
                if tmp[0] not in group_dict[tmp[1]]:
                    group_dict[tmp[1]].append(tmp[0])
        for key in group_dict:
            category_names.append(key)
            specimen_names.append(group_dict[key])
        data = {
            "task_id" : self.bind_object.sheet.id,
            "category_names": category_names,
            "specimen_names": specimen_names,
            "group_name": os.path.basename(file),
            "project_sn": self.bind_object.sheet.project_sn,
            "is_use" : 1,
        }
        col = self.db["sg_specimen_group"]
        group_id = col.insert_one(data).inserted_id
        col.update({"_id": group_id, "task_id" : self.bind_object.sheet.id}, {"$set": {"main_id": group_id}}, upsert=True)
        self.bind_object.logger.info("导入样本分组信息成功")
        return group_id, specimen_names, category_names

    def add_control_group(self,file, group_id):
        con_list = list()
        with open(file, "r") as f:
            f.readline()
            for line in f:
                tmp = line.strip().split("\t")
                if len(tmp) > 1:
                    string = tmp[0] + "|" + tmp[1]
                    con_list.append(string)
        col = self.db["sg_specimen_group"]
        result = col.find_one({"_id": group_id})
        group_name = result["group_name"]
        category_names = str(result["category_names"])
        data = {
            "task_id": self.bind_object.sheet.id,
            "compare_group_name": os.path.basename(file),
            "compare_names": json.dumps(con_list),
            "compare_category_name": "all",
            "specimen_group_id": str(group_id),
            "is_use" : 1,
        }
        col = self.db["sg_specimen_group_compare"]
        try:
            com_id = col.insert_one(SON(data)).inserted_id
            col.update({"_id": com_id, "task_id": self.bind_object.sheet.id}, {"$set": {"main_id": com_id}}, upsert=True)
        except Exception,e:
            self.bind_object.set_error("导入样本对照组信息出错:%s" % e)
        else:
            self.bind_object.logger.info("导入样本对照组信息成功")
            return com_id, con_list

if __name__ == '__main__':
    anno = Specimen(None)
    group_txt = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_itraq_and_tmt/data4/result/group.txt'
    anno.add_spe(group_txt=group_txt)
