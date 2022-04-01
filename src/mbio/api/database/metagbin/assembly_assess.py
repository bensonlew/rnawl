# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'@20181214


from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
from types import StringTypes
from biocluster.config import Config
from bson.son import SON

class AssemblyAssess(Base):
    def __init__(self, bind_object):
        super(AssemblyAssess, self).__init__(bind_object)
        self._project_type = "metagbin"

    @report_check
    def add_assembly_assess(self, genome_id, params=None, name=None):
        old_task_id = str(self.bind_object.sheet.id).split('_')[0]
        num = str(self.bind_object.sheet.id).split('_')[1]
        task_id = old_task_id + "_" + num
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "assembly_assess主表",
            "params":json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "genome_id":genome_id,# G_bin_01
            "name": name if name else "assembly_assess_Origin",
        }
        collection = self.db["assembly_assess"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_assembly_assess_detail(self, inserted_id, soft, genome_id, file_path=None):
        data_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip('\r\n').split("\t")
                if soft == "BUSCO":
                    data = {
                        "assess_id": inserted_id,
                        "genome_id": genome_id,
                        "comp": line[0],
                        "gene_dup":line[1],
                        "gene_frag": line[2],
                        "missing": line[3],
                    }
                else:
                    data = {
                        "assess_id": inserted_id,
                        "genome_id": genome_id,
                        "comp": line[0],
                        "contamina":line[1],
                        "strain_hetero": line[2],
                    }
                data_son = SON(data)
                data_list.append(data_son)
        try:
            collection = self.db["assembly_assess_detail"]
            collection.insert_many(data_list)
            main_collection = self.db["assembly_assess"]
            main_collection.update({"_id": ObjectId(inserted_id)},
                                   {"$set": {"status": 'end',
                                             'main_id':ObjectId(inserted_id),
                                            }})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (file_path, e))
        else:
            self.bind_object.logger.info("导入%s结果表成功!" % file_path)