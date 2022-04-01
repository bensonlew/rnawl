# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:20180313
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
from bson.son import SON
from bson.objectid import ObjectId
from biocluster.config import Config


class Secretory(Base):
    def __init__(self, bind_object):
        super(Secretory, self).__init__(bind_object)
        self._project_type = "fungigenome"
        self.client = Config().get_mongo_client(mtype="metagenomic", ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname("metagenomic", ref=True)]
        self.secretory = self.mongodb.secretory

    @report_check
    def add_secretory(self, out_dir, anno_summary_id=None, kegg_id=None, main=False, task_id=None, main_id=None,project_sn=None,
                      params=None, name=None, specimen_id=None):
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        nr_collection = self.db['anno_summary_detail']
        kegg_collection = self.db['anno_kegg_detail']
        data_list = []
        collection_detail = self.db['anno_secretory_detail']
        collection_type = self.db['anno_secretory_type']
        if not isinstance(anno_summary_id, ObjectId):
            anno_summary_id = ObjectId(anno_summary_id)
        if not isinstance(kegg_id, ObjectId):
            kegg_id = ObjectId(kegg_id)
        if main:
            if task_id is None:
                task_id = self.bind_object.sheet.id
            project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '分泌系统注释',
                'created_ts': created_ts,
                'name': name if name else 'Secretory_Origin',
                'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
                'status': 'end'
            }
            collection = self.db['anno_secretory']
            # 将主表名称写在这里
            main_id = collection.insert_one(insert_data).inserted_id
        else:
            if main_id is None:
                #raise Exception("main为False时需提供main_id!")
                self.bind_object.set_error("main为False时需提供main_id!", code="52102801")
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
        ko_type = {}
        type_num = {}
        kegg_type = self.secretory.find()
        for doc in kegg_type:
            ko_type[doc['KO']] = doc['type']
        with open(out_dir + '/' + specimen_id + '_secretion_system_genes.xls', 'w') as out1:
            out1.write("Gene ID\tLocation\tSample Name\tDescription\tPathway Gene Name\tKO ID\tType\n")
            for i in kegg_collection.find({"kegg_id": kegg_id, "specimen_id": specimen_id}):
                if "ko03070" in i['pathway'].split(';'):
                    data = [("anno_secretory_id", main_id)]
                    type = ko_type[i["ko"]]
                    if type in type_num:
                        type_num[type] += 1
                    else:
                        type_num[type] = 1
                    one = nr_collection.find_one(
                        {"summary_id": anno_summary_id, "specimen_id": specimen_id, "gene_id": i["gene_id"]})
                    gene_id = i["gene_id"]
                    location = i["location"]
                    gene_name = i["gene_name"]
                    ko_id = i["ko"]
                    des = one["gene_des"]
                    out1.write('\t'.join([gene_id, location, specimen_id, des, gene_name, ko_id, type]) + "\n")
                    data.extend(
                        [("gene_id", gene_id), ("specimen_id", specimen_id), ("location", location),
                         ("gene_name", gene_name), ("ko_id", ko_id), ("type", type), ("des", des)])
                    data_son = SON(data)
                    data_list.append(data_son)
        if data_list:
            try:
                collection_detail.insert_many(data_list)
                main_collection = self.db['anno_secretory']
                main_collection.update({"_id": main_id}, {"$set": {"status": 'end',"origin_id": main_id}})
            except Exception as e:
                self.bind_object.logger.error("导入secretory_detail信息出错:%s")
                self.bind_object.set_error("导入secretory_detail信息出错:%s", code="52102802")
            else:
                self.bind_object.logger.info("导入secretory_detail信息成功!")
        type_list = []
        with open(out_dir + '/' + specimen_id + '_secretion_system_type.xls', 'w') as out2:
            out2.write("Sample Name\tType\tGene No.\n")
            for j in type_num:
                data = [("anno_secretory_id", main_id)]
                number = type_num[j]
                data.extend([("specimen_id", specimen_id), ("type", j), ("anno_num", number)])
                data_son = SON(data)
                type_list.append(data_son)
                out2.write('\t'.join([specimen_id, j, str(number)]) + "\n")
        if type_list:
            try:
                collection_type.insert_many(type_list)
            except Exception as e:
                self.bind_object.logger.error("导入secretory_type信息出错:%s")
                self.bind_object.set_error("导入secretory_type信息出错:%s", code="52102803")
            else:
                self.bind_object.logger.info("导入secretory_type信息成功!")
        return main_id
