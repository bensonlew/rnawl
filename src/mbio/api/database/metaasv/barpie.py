# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.api.database.base import Base, report_check
from mainapp.libs.param_pack import group_detail_sort, param_pack
import re
import json
import datetime
from bson.objectid import ObjectId
from types import StringTypes
from mbio.packages.metaasv.common_function import find_group_name


class Barpie(Base):

    def __init__(self, bind_object):
        super(Barpie, self).__init__(bind_object)
        self._project_type = 'metaasv'
        self.task_id = ""
        self.name_id = dict()

    @report_check
    def add_barpie(self, params, from_otu_table, name=None, newick_id=None, spname_spid=None, group_id=None):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!")
        collection = self.db["asv"]
        result = collection.find_one({"_id": from_otu_table})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}barpie在表里找到相应的记录".format(str(from_otu_table)))
            self.bind_object.set_error("在sg_otu表中找不到相应记录")
        project_sn = result['project_sn']
        self.task_id = result['task_id']
        if spname_spid and params:
            if group_id not in [None, "All", "all", "ALL"]:
                ## 调用common模块，功能将导入的分组方案返回group_detail
                group_detail = find_group_name(self.task_id)
            else:
                group_detail = {'All': [str(i) for i in spname_spid.values()]}
            params['group_detail'] = group_detail_sort(group_detail)
            params['level_id'] = int(7)
        params = param_pack(params)
        if not name:
            name = "barpie_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        insert_data = {
            "project_sn": project_sn,
            'task_id': self.task_id,
            'asv_id': str(from_otu_table),
            'name': "CommunityBarPie_Origin",
            "params": params,
            'status': 'end',
            'desc': 'BarPie主表',
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        collection = self.db["barpie"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_sg_otu_detail(self, file_path, main_id):
        self.bind_object.logger.info("开始导入barpie_detail表")
        find_otu = self.db['barpie'].find_one({"_id": ObjectId(main_id)})
        if find_otu:
            self.task_id = find_otu['task_id']
        else:
            self.bind_object.set_error("没有找到相关的主表信息")
        spe_str=""  #guanqing.zou 物种按丰度排序，以|分割的字符串 20180411
        with open(file_path, 'rb') as r:
            lines = r.readlines()
            head = lines[0].strip('\r\n')
            head = re.split('\t', head)
            sample_list = head[1:]
            data_list = []
            rank = 1
            for line in lines[1:]:
                line = line.strip().split('\t')
                sample_num = [float(x) for x in line[1:]]
                classify_list = re.split(r"\s*;\s*", line[0])
                spe_str+=classify_list[-1]+'|'       #guanqing 20180411
                if len(re.split("; ", line[0])) > 1:
                    species_name = line[0].strip().split("; ")[-1].strip()
                elif len(re.split(";", line[0])) > 1:
                    species_name = line[0].strip().split(";")[-1].strip()
                else:
                    species_name = line[0]
                otu_detail = {
                    "barpie_id": ObjectId(main_id),
                    "species_name": species_name,
                    "rank": rank
                }
                rank += 1
                values_dict = dict(zip(sample_list, sample_num))
                data_list.append(dict(otu_detail, **values_dict))
            try:
                collection = self.db['barpie_detail']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.logger.error("导入barpie_detail表格信息出错:{}".format(e))
            else:
                self.bind_object.logger.info("导入barpie_detail表格成功")
            try:
                main_collection = self.db['barpie']
                settled_params = {'software': "python-2.7"}
                settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
                column_data = {
                    "column_data": {"name":"species_name",
                            "data": sample_list,
                            "category": ""}
                    }
                column_data_json = json.dumps(column_data, sort_keys=True, separators=(',', ':'))
                pie_data = {
                    "pie_data": {"name":"species_name",
                            "data": sample_list,
                            "category": ""}
                    }
                pie_data_json = json.dumps(pie_data, sort_keys=True, separators=(',', ':'))
                main_collection.update({"_id":ObjectId(main_id)},{"$set":{"settled_params":settled_params_json,
                                                                          "main_id": ObjectId(main_id),
                                                                          "column_data":column_data_json,
                                                                          "pie_data":pie_data_json,
                                                                          }})
            except Exception as e:
                self.bind_object.logger.error("更新barpie表格信息出错:{}".format(e))
            else:
                self.bind_object.logger.info("更新barpie表格成功")
