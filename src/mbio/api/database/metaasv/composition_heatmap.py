# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.api.database.base import Base, report_check
from mainapp.libs.param_pack import group_detail_sort, param_pack
import json
import datetime
from bson.objectid import ObjectId
from types import StringTypes
import re
from mbio.packages.metaasv.common_function import find_group_name


class CompositionHeatmap(Base):
    def __init__(self, bind_object):
        super(CompositionHeatmap, self).__init__(bind_object) #
        self._project_type = 'metaasv'
        self.task_id = ""
        self.name_id = dict()
        self.main_task_id = "_".join(self.bind_object.sheet.id.split('_')[0:2])  # add task_id by guhaidong 20171120

    @report_check
    def add_heatmap(self, params, from_otu_table, name=None, spname_spid=None,group_id=None):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!")
        collection = self.db["asv"]
        result = collection.find_one({"_id": from_otu_table})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}在asv表里找到相应的记录".format(str(from_otu_table)))
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
            name = "CommunityHeatmap_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S") #
        insert_data = {
            "project_sn": project_sn,
            'task_id': self.task_id,
            'asv_id': from_otu_table,
            'name': "CommunityHeatmap_Origin",
            "params": params,
            'status': 'end',
            'desc': 'Composition_Heatmap主表', #
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        collection = self.db["heatmap"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_heatmap_detail(self, file_path, main_id,type, species_sorts=None, specimen_sorts=None):
        if main_id != 0 and not isinstance(main_id, ObjectId):
            if isinstance(main_id, StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error("main_id必须为ObjectId对象或其对应的字符串!")
        self.bind_object.logger.info("开始导入heatmap_detail表")
        data_list = []
        species_sort = []
        with open(file_path, 'rb') as r:
            head = r.next().strip('').split('\t')
            specimen_list = head[1:-1]
            index = 0
            for line in r:
                line = line.strip().split('\t')
                sample_num = line[1:-1]
                index += 1
                if len(re.split("; ", line[0])) > 1:
                    species_name = line[0].strip().split("; ")[-1].strip()
                elif len(re.split(";", line[0])) > 1:
                    species_name = line[0].strip().split(";")[-1].strip()
                else:
                    species_name = line[0]
                if species_name not in species_sort:
                    species_sort.append(species_name)
                otu_detail = {
                    "heatmap_id": ObjectId(main_id),
                    "species_name": species_name,
                    "type": type,
                    "level_color": line[-1]
                    }
                otu_detail['rank']=index
                values_dict = dict(zip(specimen_list, sample_num))
                data_list.append(dict(otu_detail, **values_dict))
        try:
            collection = self.db['heatmap_detail']
            collection.insert_many(data_list)
            self.bind_object.logger.info("已写入")
        except Exception as e:
            self.bind_object.logger.error("导入heatmap_detail表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入heatmap_detail表格成功")
        try:
            main_collection = self.db['heatmap']
            settled_params = {'software': "R-3.3.1, python-2.7"}
            settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
            heatmap_data = {
                "heatmap_data": {"name":"species_name",
                        "data": specimen_list}
                }
            heatmap_data_json = json.dumps(heatmap_data, sort_keys=True, separators=(',', ':'))
            tree_data = {
                "tree_data": {"name":"name",
                        "condition": {"type":["species", "specimen"]}}
                }
            tree_data_json = json.dumps(tree_data, sort_keys=True, separators=(',', ':'))
            if species_sorts:
                species_sort = species_sorts
            else:
                species_sort = species_sort
            if specimen_sorts:
                specimen_sort = specimen_sorts
            else:
                specimen_sort = specimen_list
            main_collection.update_one({"_id": ObjectId(main_id)},{"$set": {"main_id": ObjectId(main_id),
                                                                "species_sort": species_sort,
                                                                "specimen_sort": specimen_sort,
                                                                "settled_params": settled_params_json,
                                                                "heatmap_data": heatmap_data_json,
                                                                "tree_data": tree_data_json,}})
        except Exception as e:
            self.bind_object.logger.error("更新heatmap表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("更新heatmap表格成功")

    def insert_tree_table(self, file_path, update_id, type):
        if update_id is None:
            self.bind_object.set_error("需提供dist_id!")
        else:
            if not isinstance(update_id, ObjectId):
                dist_id = ObjectId(update_id)
            else:
                dist_id = update_id
        data_list = []
        with open(file_path, 'r') as f:
            line = f.readline().strip()
            insert_data = {
                "name": "",
                "data": line,
                "type": type,
                "heatmap_id":dist_id
                }
            data_list.append(insert_data)
        try:
            collection = self.db["heatmap_tree"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("更新%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("更新%s信息成功!" % file_path)