# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from biocluster.api.database.base import Base, report_check
import datetime
import json,os,re
from bson.son import SON
import pandas as pd
from bson.objectid import ObjectId

class CoreGene(Base):
    def __init__(self, bind_object):
        super(CoreGene, self).__init__(bind_object)
        self._project_type = "metagbin"

    @report_check
    def add_coregene(self, params=None, name=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "Core_gene",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": name if name else "Core_gene"
        }
        collection = self.db["core_gene"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_coregene_detail(self, inserted_id, dir):
        files = os.listdir(dir)
        for file in files:
            if re.search(r'.xls', file):
                data_list = []
                with open(dir + '/' + file, 'r') as f:
                    lines = f.readlines()
                    for lin in lines[1:]:
                        line = lin.rstrip('\r\n').split('\t')
                        data = {
                            "cor_id": ObjectId(inserted_id),
                            "bin_id": line[0],
                            "location": line[1],
                            "start": line[2],
                            "end": line[3],
                            "coregene_name": line[4],
                            "indentity": line[5],
                            "coverage": line[6],
                            "seq": line[7],
                        }
                        data_son = SON(data)
                        data_list.append(data_son)
                try:
                    collection = self.db["core_gene_detail"]
                    collection.insert_many(data_list)
                except Exception, e:
                    self.bind_object.logger.info("导入%s结果表出错:%s" % (file, e))
                else:
                    self.bind_object.logger.info("导入%s结果表成功!" % file)
            elif re.search(r'.fa', file):
                data2_list = []
                bin_id = file.split('.')[0]
                with open(dir + '/' + file, 'r') as f:
                    lines = f.readlines()
                    seq = lines[1].rstrip('\r\n')
                    data ={
                        "cor_id": ObjectId(inserted_id),
                        "bin_id": bin_id,
                        "seq": seq
                    }
                    data_son = SON(data)
                    data2_list.append(data_son)
                try:
                    collection = self.db["core_gene_seq"]
                    collection.insert_many(data2_list)
                except Exception, e:
                    self.bind_object.logger.info("导入%s结果表出错:%s" % (file, e))
                else:
                    self.bind_object.logger.info("导入%s结果表成功!" % file)