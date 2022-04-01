# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.api.database.base import Base, report_check
import re
import json
import datetime
from bson.objectid import ObjectId
from types import StringTypes


class Circos(Base):

    def __init__(self, bind_object):
        super(Circos, self).__init__(bind_object)
        self._project_type = 'metaasv'
        self.task_id = ""
        self.name_id = dict()

    @report_check
    def add_barpie(self, params, from_otu_table, name=None, newick_id=None):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!")
        collection = self.db["asv"]
        result = collection.find_one({"_id": from_otu_table})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}circos在表里找到相应的记录".format(str(from_otu_table)))
            self.bind_object.set_error("在sg_otu表中找不到相应记录")
        project_sn = result['project_sn']
        self.task_id = result['task_id']
        if not name:
            name = "circos_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        insert_data = {
            "project_sn": project_sn,
            'task_id': self.task_id,
            'asv_id': str(from_otu_table),
            'name': self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else name,
            "params": params,
            'status': 'end',
            'desc': 'Circos主表',
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        collection = self.db["circos"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_sg_otu_detail(self, file_path, main_id):
        self.bind_object.logger.info("开始导入circos_detail表")
        find_otu = self.db['circos'].find_one({"_id": ObjectId(main_id)})
        if find_otu:
            self.task_id = find_otu['task_id']
        else:
            self.bind_object.set_error("没有找到相关的主表信息")
        species_sort = []
        with open(file_path, 'rb') as r:
            lines = r.readlines()
            head = lines[0].strip('\r\n')
            head = re.split('\t', head)
            sample_list = head[1:]
            data_list = []
            for line in lines[1:]:
                line = line.strip().split('\t')
                sample_num = [float(x) for x in line[1:]]
                if line[0] not in species_sort:
                    species_sort.append(line[0].strip().split("; ")[-1].strip())
                otu_detail = {
                    "circos_id": ObjectId(main_id),
                    "species_name": line[0].strip().split("; ")[-1].strip()
                }
                values_dict = dict(zip(sample_list, sample_num))
                data_list.append(dict(otu_detail, **values_dict))
        try:
            collection = self.db['circos_detail']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入circos_detail表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入circos_detail表格成功")
        try:
            main_collection = self.db['circos']
            settled_params = {'software': "python-2.7"}
            settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
            specimen_sort = sample_list
            main_collection.update({"_id":ObjectId(main_id)},{"$set":{"settled_params":settled_params_json,
                                                                      "specimen_sort":specimen_sort,
                                                                      "main_id": ObjectId(main_id),
                                                                      "species_sort":species_sort,
                                                                      }})
        except Exception as e:
            self.bind_object.logger.error("更新circos表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("更新circos表格成功")
