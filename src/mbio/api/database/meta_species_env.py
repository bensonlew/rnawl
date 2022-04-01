# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
# from biocluster.config import Config
from mainapp.libs.param_pack import group_detail_sort
import os

class MetaSpeciesEnv(Base):
    def __init__(self, bind_object):
        super(MetaSpeciesEnv, self).__init__(bind_object)
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB

    @report_check
    def add_mantel_table(self, level, otu_id, env_id, task_id=None, name=None, params=None, spname_spid=None):
        if level not in range(1, 10):
            self.bind_object.logger.error("level参数%s为不在允许范围内!" % level)
            self.bind_object.set_error("level参数不在范围内", code="51003901")
        # if task_id is None:
        #     task_id = self.bind_object.sheet.id
        if name:
            origin_name = name
        else:
            origin_name = "mantel_origin"
        matrix_types = ["species_matrix", "env_matrix"]
        if "units" in params:
            matrix_types = ["species_matrix", "env_matrix", "partial_matrix"]
            origin_name = "partial_mantel_origin"

        if not isinstance(otu_id, ObjectId):
            otu_id = ObjectId(otu_id)
        collection = self.db["sg_otu"]
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
            "otu_id": otu_id,
            "name": self.bind_object.sheet.main_table_name if self.bind_object.sheet.main_table_name else origin_name,
            "level_id": level,
            "status": "end",
            "desc": "",
            "matrix_types": matrix_types,
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_species_mantel_check"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_mantel_detail(self, file_path, mantel_id=None, main_colletion_name=None):
        main_colletion_time = ""
        if main_colletion_name:
            main_colletion_time = main_colletion_name[11:]
        data_list = []
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
                        data["permutations"] = line[6]
                        data["mantel_r"] = line[4]
                        data["p_value"] = line[5]
                        data["tail_type"] = line[7]
                        data["cdm"] = "partial_matrix"+main_colletion_time
                    else:
                        data["entries_num"] = line[2]
                        data["permutations"] = line[5]
                        data["mantel_r"] = line[3]
                        data["p_value"] = line[4]
                        data["tail_type"] = line[6]
                    data_list.append(data)
        try:
            collection = self.db["sg_species_mantel_check_detail"]
            collection.insert_many(data_list)

            collection_main = self.db["sg_species_mantel_check"]
            #collection_main.update({"_id": mantel_id}, {"$set": {"main_id" : mantel_id}})
        except Exception, e:
            self.bind_object.logger.error("导入mantel检验结果数据出错:%s" % e)
            self.bind_object.set_error("导入mantel检验结果数据出错", code="51003902")
        else:
            self.bind_object.logger.info("导入mantel检验结果数据成功")

    @report_check
    def add_mantel_matrix(self, matrix_path, matrix_type, mantel_id=None):
        data_list = []
        with open(matrix_path, "r") as m:
            samples = m.readline().strip().split("\t")
            # sample_num = len(samples)
            # self.bind_object.logger.info(samples)
            for line in m:
                line = line.strip().split("\t")
                data = {
                    "mantel_id": mantel_id,
                    "matrix_type": matrix_type,
                    "column_name": line[0]
                }
                for n, s in enumerate(samples):
                    data[s] = line[n+1]
                data_list.append(data)
        # self.bind_object.logger.info data_list
        try:
            collection = self.db["sg_species_mantel_check_matrix"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入mantel检验%s矩阵出错:%s" % matrix_type,e)
            self.bind_object.set_error("导入mantel检验矩阵出错", code="51003903")
        else:
            self.bind_object.logger.info("导入mantel检验%s矩阵成功" % matrix_type)

    @report_check
    def add_correlation(self, level, otu_id, env_id, species_tree=None, env_tree=None, task_id=None, name=None, params=None, spname_spid=None, env_list=None, species_list=None):
        if level not in range(1, 10):
            self.bind_object.logger.error("level参数%s为不在允许范围内!" % level)
            self.bind_object.set_error("level参数不在范围内", code="51003901")
        # if task_id is None:
        #     task_id = self.bind_object.sheet.id
        if not isinstance(otu_id, ObjectId):
            otu_id = ObjectId(otu_id)
        collection = self.db["sg_otu"]
        result = collection.find_one({"_id": otu_id})
        if task_id is None:
            task_id = result['task_id']
        if spname_spid:
            group_detail = {'All': [str(i) for i in spname_spid.values()]}
            params['group_detail'] = group_detail_sort(group_detail)
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": task_id,
            "env_id": ObjectId(env_id),
            "otu_id": otu_id,
            "name": name if name else "pearson_origin",
            "level_id": level,
            "status": "end",
            "env_tree": env_tree if env_tree else "()",
            "species_tree": species_tree if species_tree else "()",
            "env_list": env_list if env_list else "[]",
            "species_list": species_list if species_list else "[]",
            "desc": "",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["sg_species_env_correlation"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_correlation_detail(self, file_path, value_type, correlation_id=None, species_tree=None, env_tree=None, env_list=None, species_list=None,sorted_otu_file=None):
        data_list = []

        ###guanqing.zou 20180419 生成分类水平的丰度值按高到底排序的列表

        if sorted_otu_file:
            sort_otu_path = sorted_otu_file
        else:
            sort_otu_path='/'.join(file_path.split('/')[:-2])+'/sorted_otu_file.xls'
        if not os.path.exists(sort_otu_path):
            self.bind_object.logger.error(sort_otu_path+' 文件不存在')
            self.bind_object.set_error("sort_otu_path文件不存在", code="51003904")
        sort_spe_list=[]
        with open(sort_otu_path,'r') as f1:
            for line in f1:
                sort_spe_list.append(line.split('\t')[0].split(';')[-1].strip())
        #####Done


        with open(file_path, "r") as f:
            envs = f.readline().strip().split("\t")[1:]
            for line in f:
                line = line.strip().split("\t")
                data = {
                    "correlation_id": ObjectId(correlation_id),
                    "species_name": line[0],
                    "value_type": value_type,
                    "rank": sort_spe_list.index(line[0])   #guanqing.zou 20180419
                }
                for n, e in enumerate(envs):
                    data[e] = line[n+1]
                data_list.append(data)
        try:
            collection = self.db["sg_species_env_correlation_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["sg_species_env_correlation"]
            main_collection.update({"_id": ObjectId(correlation_id)},
                                   {"$set": {
                                    "env_tree": env_tree if env_tree else "()",
                                    "species_tree": species_tree if species_tree else "()",
                                    "env_list": env_list if env_list else "[]",
                                    "is_sort" : 1,
                                    "species_list": species_list if species_list else "[]"}})
            #main_collection.update({"_id": ObjectId(correlation_id)}, {"$set": {"main_id" :  ObjectId(correlation_id)}})
        except Exception, e:
            self.bind_object.logger.error("导入皮尔森相关系数矩阵出错:%s" % e)
            self.bind_object.set_error("导入皮尔森相关系数矩阵出错", code="51003905")
        else:
            self.bind_object.logger.info("导入皮尔森相关系数矩阵成功")
