# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json


class Vpa(Base):
    def __init__(self, bind_object):
        super(Vpa, self).__init__(bind_object)
        self._project_type = "metaasv"

    @report_check
    def add_vpa(self,level_id,otu_id, name=None, params=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "level_id": level_id,
            "otu_id" : otu_id,
            "status": "end",
            "desc": "正在计算",
            "name": name,
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["vpa"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_vpa_graph(self, vpa_graph_path, vpa_id):
        data_list = []
        if not isinstance(vpa_id, ObjectId):
            vpa_id = ObjectId(vpa_id)
        with open(vpa_graph_path, "r") as f:
            f.readline()
            for line in f:
                line = line.strip().split("\t")
                self.bind_object.logger.info("line[0]: %s\n" % line[0])
                if line[0].strip() in ["Residuals"]:
                    insert_data = {
                        'vpa_id': vpa_id,
                        'type' : 'residual'
                        }
                    try:
                        insert_data['num'] = float(line[1])
                    except:
                        insert_data['num'] = line[1]##结果为NA
                    data_list.append(insert_data)
                else:
                    insert_data = {
                        'vpa_id': vpa_id,
                        'label': line[0],
                        'type' : 'r2'
                        }
                    try:
                        insert_data['num'] = float(line[1])
                    except:
                        insert_data['num'] = line[1]##结果为NA
                    data_list.append(insert_data)
        try:
            collection = self.db["vpa_venn"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.error("导入Vpa画图数据出错:%s" % e)
            self.bind_object.set_error("import Vpa graph data ERROR")
        else:
            self.bind_object.logger.info("导入Vpa画图数据成功")
        try:
            main_collection = self.db["vpa"]
            settled_params = {"software" : "R-3.3.1 (vegan)"}
            settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
            venn_data = {"venn_data": {"names":"num","category":"label"}}
            venn_data_json = json.dumps(venn_data, sort_keys=True, separators=(',', ':'))
            text_data = {"text_data": {"name":"num", "condition":{"type": "residual"}}}
            text_data_json = json.dumps(text_data, sort_keys=True, separators=(',', ':'))
            main_collection.update({"_id": ObjectId(vpa_id)},
                                   {"$set": {"status": 'end',
                                            "main_id": vpa_id,
                                             "settled_params": settled_params_json,
                                             "venn_data": venn_data_json,
                                             "text_data": text_data_json}})
        except Exception, e:
            self.bind_object.logger.error("更新Vpa画图数据出错:%s" % e)
            self.bind_object.set_error("import Vpa graph data ERROR")
        else:
            self.bind_object.logger.info("更新Vpa画图数据成功")

    @report_check
    def add_vpa_detail(self, vpa_path, vpa_id):
        data_list = []
        if not isinstance(vpa_id, ObjectId):
            vpa_id = ObjectId(vpa_id)
        with open(vpa_path, 'rb') as r:
                r.readline()
                for line in r:
                    line = line.strip().split("\t")
                    insert_data = {
                        'vpa_id': vpa_id,
                        'env' : line[0],
                        'radj': line[1],
                        'group' : line[2]
                        }
                    data_list.append(insert_data)
        try:
            collection = self.db["vpa_table"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入Vpa数据出错:%s" % e)
            self.bind_object.set_error("import Vpa data ERROR")
        else:
            self.bind_object.logger.info("导入Vpa数据成功")
        try:
            main_collection = self.db["vpa"]
            table_data = {"table_data": ["env", "group", "radj"]}
            table_data_json = json.dumps(table_data, sort_keys=True, separators=(',', ':'))
            main_collection.update({"_id":vpa_id},{"$set": {"main_id": vpa_id,
                                                            "table_data": table_data_json}})
        except Exception, e:
            self.bind_object.logger.info("更新Vpa数据出错:%s" % e)
            self.bind_object.set_error("import Vpa data ERROR")
        else:
            self.bind_object.logger.info("更新Vpa数据成功")
