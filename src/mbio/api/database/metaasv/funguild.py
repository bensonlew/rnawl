# -*- coding: utf-8 -*-
# __author__ = "qingchen.zhang"

from biocluster.api.database.base import Base, report_check
import json
import datetime
from bson.objectid import ObjectId
import pandas as pd


class Funguild(Base):
    def __init__(self, bind_object):
        super(Funguild, self).__init__(bind_object)
        self._project_type = "metaasv"

    @report_check
    def add_funguild(self,level_id,otu_id, name=None, params=None):
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
        collection = self.db["funguild"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_detail(self, file_path,table_id):
        data_list = []
        if not isinstance(table_id, ObjectId):
            table_id = ObjectId(table_id)
        data = pd.read_table(file_path, sep="\t")
        data.fillna("-",inplace=True)
        taxon_list = ["taxon_" + i for i in "d,k,p,c,o,f,g,s".split(",")]
        head_list = data.columns.tolist()
        sample_list = head_list[1:head_list.index("taxonomy")]
        for i in range(len(data)):
            insert_data = {
                "guild_id": table_id,
                "asv_name": data["ASV ID"][i],
                "guild": data["Guild"][i],
                "confidence" : data["Confidence Ranking"][i],
                "growth" : data["Growth Morphology"][i],
                "trait" : data["Trait"][i] ,
                "note" : data["Notes"][i],
                "source" : data["Citation/Source"][i],
                "trophic":data["Trophic Mode"][i]
            }
            for sample in sample_list:
                insert_data[sample] = data[sample][i]

            split_taxon = data["taxonomy"][i].split(";")
            for id,taxon in enumerate(taxon_list,0):
                try:
                    insert_data[taxon]= split_taxon[id].strip()
                except Exception as e:
                    insert_data[taxon]= "-"

            data_list.append(insert_data)
        try:
            collection = self.db["funguild_detail"]
            collection.insert_many(data_list)
        except Exception as  e:
            self.bind_object.logger.info("导入Funguild 详情数据出错:%s" % e)

        else:
            self.bind_object.logger.info("导入Funguild 详情数据成功")
        try:
            main_collection = self.db["funguild"]
            settled_params = {"software" : "FUNGuild v1.0"}
            settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(",", ":"))
            table_data = {"table_data": ["asv_name"]+sample_list + taxon_list + ["guild","confidence","growth","trait","note","source","trophic"]}
            table_data_json = json.dumps(table_data, sort_keys=True, separators=(",", ":"))
            main_collection.update({"_id": ObjectId(table_id)},{"$set": {"status": "end",
                                                                         "main_id": ObjectId(table_id),
                                                                         "settled_params": settled_params_json,
                                                                         "table_data":table_data_json}})
        except Exception as  e:
            self.bind_object.logger.info("导入Funguild 详情数据出错:%s" % e)

        else:
            self.bind_object.logger.info("导入Funguild 详情数据成功")

    @report_check
    def add_sum(self,file_path, table_id):
        self.bind_object.logger.info("table_id:{}".format(table_id))
        if not isinstance(table_id, ObjectId):
            new_table_id = ObjectId(table_id)
        else:
            new_table_id = table_id
        data = pd.read_table(file_path,sep="\t",header=0)
        data.fillna("-",inplace=True)
        samples = data.columns[1:]
        insert_data = []
        for index in range(len(data)):
            tmp_data = {
                "guild_id": new_table_id,
                "guild" : data["Guild"][index]
            }
            for sample in samples:
                tmp_data[sample] = data[sample][index]
            insert_data.append(tmp_data)
        try:
            collection = self.db["funguild_bar"]
            collection.insert_many(insert_data)
            self.bind_object.logger.info("导入Funguild_bar统计数据成功")
        except Exception as e:
            self.bind_object.logger.error("导入funguild_bar统计数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入Funguild统计数据成功")
        try:
            main_collection = self.db["funguild"]
            column_data = {"column_data": {"name":"guild","data": list(samples),"category": ""}}
            column_data_json = json.dumps(column_data, sort_keys=True, separators=(",", ":"))
            guild_table_data = {"table_data":["guild"]+list(samples), "condition":{"type": "guild"}}
            table_data_json = json.dumps(guild_table_data, sort_keys=True, separators=(",", ":"))
            main_collection.update({"_id": new_table_id},{"$set":{"specimen_list": ",".join(list(samples)),
                                                                        "column_data": column_data_json,
                                                                        "guild_table_data": table_data_json}})
        except Exception as e:
            self.bind_object.logger.error("更新Funguild统计数据出错:%s" % e)

