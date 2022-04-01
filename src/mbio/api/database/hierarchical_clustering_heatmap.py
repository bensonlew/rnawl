# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
# last modified by guhaidong 20171120

from biocluster.api.database.base import Base, report_check
import re
import datetime
from bson.objectid import ObjectId
from types import StringTypes
# from biocluster.config import Config


class HierarchicalClusteringHeatmap(Base):
    """
    聚类不聚类heatmap图合并api接口（因为分组求的问题这边的表不能再被利用）
    """
    def __init__(self, bind_object):
        super(HierarchicalClusteringHeatmap, self).__init__(bind_object) #
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB
        self.task_id = ""
        self.name_id = dict()
        self.main_task_id = "_".join(self.bind_object.sheet.id.split('_')[0:2])  # add task_id by guhaidong 20171120

    @report_check
    def add_sg_hc_heatmap(self, params, from_otu_table, name=None, sample_tree=None, sample_list=None, species_tree=None, species_list=None):
        if from_otu_table != 0 and not isinstance(from_otu_table, ObjectId):
            if isinstance(from_otu_table, StringTypes):
                from_otu_table = ObjectId(from_otu_table)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!", code="51003601")
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": from_otu_table})
        if not result:
            self.bind_object.logger.error("无法根据传入的_id:{}在sg_otu表里找到相应的记录".format(str(from_otu_table)))
            self.bind_object.set_error("在sg_otu表中找不到相应记录", code="51003602")
        project_sn = result['project_sn']
        self.task_id = result['task_id']
        if not name:
            name = "hierarchical_clustering_heatmap_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S") #
        insert_data = {
            "project_sn": project_sn,
            'task_id': self.task_id,
            'otu_id': from_otu_table,
            'name': self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else name,
            "params": params,
            "sample_tree": sample_tree,
            "sample_list": sample_list,
            "species_tree": species_tree,
            "species_list": species_list,
            'status': 'end',
            'desc': 'otu table after Hierarchical Clustering Heatmap', #
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "show": 0,
            "type": "otu_HierarchicalClusteringHeatmap", #
            "is_sort" : 1
        }
        collection = self.db["sg_hc_heatmap"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_sg_hc_heatmap_detail(self, file_path, color_dict, new_otu_id, from_otu_id, sample_tree=None,
                                 sample_list=None, species_tree=None, species_list=None,
                                 otu_relative=None):  # add color_dict by houshuang 20190918
        if from_otu_id != 0 and not isinstance(from_otu_id, ObjectId):
            if isinstance(from_otu_id, StringTypes):
                from_otu_id = ObjectId(from_otu_id)
            else:
                self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!", code="51003603")
        self.bind_object.logger.info("开始导入sg_hc_heatmap_detail表")
        self.color_dict, self.new_otu_id = color_dict, new_otu_id
        insert_data = list()
        self.get_insertdate(insert_data, file_path)
        self.get_insertdate(insert_data, otu_relative, "relative")
        try:
            collection = self.db['sg_hc_heatmap_detail']
            collection.insert_many(insert_data)
            main_collection = self.db["sg_hc_heatmap"]
            self.bind_object.logger.info("开始刷新主表写树")
            self.bind_object.logger.info("_id is : %s; task_id is : %s" %(new_otu_id, self.main_task_id))
            main_collection.update({"_id": ObjectId(new_otu_id), "task_id": self.main_task_id},
                                    {"$set": {
                                        "sample_tree": sample_tree if sample_tree else "()",
                                        "species_tree": species_tree if species_tree else "()",
                                        "sample_list": sample_list if sample_list else [],
                                        "species_list": species_list if species_list else []}})  #  add task_id by guhaidong 20171115

            #main_collection.update({"_id":ObjectId(new_otu_id)},{"$set": {"main_id":ObjectId(new_otu_id)}})

            self.bind_object.logger.info("已写入")
        except Exception as e:
            self.bind_object.logger.error("导入sg_hc_heatmap_detail表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入sg_hc_heatmap_detail表格成功")

    def get_insertdate(self, insert_data, infile, data_type='absolute'):
        with open(infile, 'rb') as r:
            head = r.next().strip('\r\n')   #windows换行符
            head = re.split('\t', head)
            new_head = head[1:]
            index=0  # by guanqing.zou
            for line in r:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                sample_num = line[1:]
                classify_list = re.split(r"\s*;\s*", line[0])
                otu_detail = dict()
                index+=1   # by guanqing.zou 20180403
                otu_detail['rank']=index    # by guanqing.zou 20180403
                otu_detail['hc_id'] = ObjectId(self.new_otu_id)
                for cf in classify_list:
                    if cf != "":
                        otu_detail[cf[0:3].lower()] = cf
                        if self.color_dict != {}:
                            otu_detail["level_color"] = self.color_dict[cf]  # by houshuang 20191010
                for i in range(0, len(sample_num)):
                    otu_detail[new_head[i]] = sample_num[i]
                otu_detail['type'] = data_type
                # otu_detail['task_id'] = self.task_id
                insert_data.append(otu_detail)
        