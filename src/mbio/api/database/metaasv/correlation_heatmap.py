# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
from mainapp.libs.param_pack import group_detail_sort
import os

class CorrelationHeatmap(Base):
    def __init__(self, bind_object):
        super(CorrelationHeatmap, self).__init__(bind_object)
        self._project_type = 'metaasv'

    @report_check
    def add_correlation(self, level, otu_id, env_id, species_tree=None, env_tree=None, task_id=None, name=None, params=None, spname_spid=None, env_list=None, species_list=None):
        if level not in range(1, 10):
            self.bind_object.logger.error("level参数%s为不在允许范围内!" % level)
            self.bind_object.set_error("level参数不在范围内")
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
            "env_id": ObjectId(env_id),
            "asv_id": otu_id,
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
        collection = self.db["cor_heatmap"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_correlation_detail(self, file_path, value_type, correlation_id=None, species_tree=None, env_tree=None, env_list=None, species_list=None):
        data_list = []

        ###guanqing.zou 20180419 生成分类水平的丰度值按高到底排序的列表
        
        sort_otu_path='/'.join(file_path.split('/')[:-2])+'/abundance.xls'
        if not os.path.exists(sort_otu_path):
            self.bind_object.logger.error(sort_otu_path+' 文件不存在')
            self.bind_object.set_error("sort_otu_path文件不存在")
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
                    "heatmap_id": ObjectId(correlation_id),
                    "species_name": line[0],
                    "type": value_type,
                    "rank": sort_spe_list.index(line[0])   #guanqing.zou 20180419
                }
                for n, e in enumerate(envs):
                    try:
                        data[e] = float(line[n+1])
                    except:
                        data[e] = line[n+1]
                data_list.append(data)
        try:
            collection = self.db["cor_heatmap_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入皮尔森相关系数矩阵出错:%s" % e)
            self.bind_object.set_error("导入皮尔森相关系数矩阵出错")
        else:
            self.bind_object.logger.info("导入皮尔森相关系数矩阵成功")
        try:
            main_collection = self.db['cor_heatmap']
            settled_params = {'software': "R-3.3.1, python-2.7"}
            settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
            heatmap_data = {
                "heatmap_data": {"name":"species_name",
                    "data": envs}
                }
            heatmap_data_json = json.dumps(heatmap_data, sort_keys=True, separators=(',', ':'))
            table_data = {
                "table_data": {"name":"species_name",
                    "data": envs, "condition": {"type": ["correlation", "pvalue"]}}
                }
            table_data_json = json.dumps(table_data, sort_keys=True, separators=(',', ':'))
            tree_data = {
                "tree_data": {"name":"name",
                        "condition": {"type":["species", "specimen"]}}
                }
            tree_data_json = json.dumps(tree_data, sort_keys=True, separators=(',', ':'))
            if env_list:
                species_sort = env_list
            else:
                species_sort = []
            if species_list:
                specimen_sort = species_list
            else:
                specimen_sort = []
            main_collection.update_one({"_id": ObjectId(correlation_id)},{"$set": {"main_id": ObjectId(correlation_id),
                                                                "species_sort": species_sort,
                                                                "specimen_sort": specimen_sort,
                                                                "table_data": table_data_json,
                                                                "settled_params": settled_params_json,
                                                                "heatmap_data": heatmap_data_json,
                                                                "tree_data": tree_data_json}})
        except Exception as e:
            self.bind_object.logger.error("更新cor_heatmap表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("更新cor_heatmap表格成功")

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
            collection = self.db["cor_heatmap_tree"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("更新%s信息出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("更新%s信息成功!" % file_path)