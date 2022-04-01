# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# last modify: 2017.10.17

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
from types import StringTypes
# from biocluster.config import Config
from mainapp.libs.param_pack import group_detail_sort


class HeatmapCor(Base):
    def __init__(self, bind_object):
        super(HeatmapCor, self).__init__(bind_object)
        self._project_type = "metagenomic"
        # self._db_name = Config().MONGODB + '_metagenomic'

    @report_check
    def add_heatmap_cor(self, anno_type, name=None, params=None):
        task_id = self.bind_object.sheet.id
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": task_id,
            "anno_type": anno_type,
            #"env_id": ObjectId(env_id),
            #"geneset_id": ObjectId(geneset_id),
            #"anno_id": anno_id,
            "name": name,
            #"level_id": level_id,
            "status": "end",
            #'group_id': group_id,
            "desc": "相关性热图",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["heatmap_cor"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_heatmap_cor_detail(self, file_path, value_type, heatmap_cor_id, species_tree=None, env_tree=None):
        env_list = ""
        species_list = ""
        if env_tree != "":
            with open(env_tree, 'r') as f:
                env_tree = f.readline()
        if species_tree != "":
            with open(species_tree, 'r') as f:
                species_tree = f.readline()
        data_list = []
        with open(file_path, "r") as f:
            envs = f.readline().strip().split("\t")[1:]
            env_list = ','.join(envs)
            for line in f:
                line = line.strip().split("\t")
                data = {
                    "heatmap_cor_id": ObjectId(heatmap_cor_id),
                    "species_name": line[0],
                    "type": value_type
                }
                species_list = species_list + ";" + line[0]
                for n, e in enumerate(envs):
                    data[e] = line[n + 1]
                data_list.append(data)
            species_list = species_list.lstrip(';')
        try:
            collection = self.db["heatmap_cor_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["heatmap_cor"]
            main_collection.update({"_id": ObjectId(heatmap_cor_id)},
                                   {"$set": {
                                       "env_tree": env_tree,
                                       "species_tree": species_tree,
                                       "env_list": env_list,
                                       "species_list": species_list}})
        except Exception, e:
            self.bind_object.logger.info("导入皮尔森相关系数矩阵出错:%s" % e)
            self.bind_object.set_error("导入皮尔森相关系数矩阵出错", code="52801301")
        else:
            self.bind_object.logger.info("导入皮尔森相关系数矩阵成功")
