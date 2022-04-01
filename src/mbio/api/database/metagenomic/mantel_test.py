# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'


from biocluster.api.database.base import Base, report_check
import json
import datetime
from bson.objectid import ObjectId
#from biocluster.config import Config


class MantelTest(Base):
    def __init__(self, bind_object):
        super(MantelTest, self).__init__(bind_object)
        self._project_type = "metagenomic"
        #self._db_name = Config().MONGODB + '_metagenomic'

    @report_check
    def add_mantel_test(self, anno_type, name=None, params=None):
        task_id = self.bind_object.sheet.id
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": task_id,
            "anno_type": anno_type,
            #"env_id": env_id,
            #"geneset_id": geneset_id,
            #"anno_id": anno_id,
            "name": name,
            #"level_id": level_id,
            "status": "end",
            #'group_id': group_id,
            "desc": "mantel test 分析",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        collection = self.db["mantel_test"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_mantel_test_detail(self, file_path, mantel_id):
        data_list = []
        with open(file_path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                elif line.startswith("DM1"):
                    continue
                else:
                    line = line.strip().split("\t")
                    data = {
                        "mantel_id": ObjectId(mantel_id),
                        "dm1": "species_matrix",
                        "dm2": "env_matrix"
                    }
                    if len(line) == 8:
                        data["entries_num"] = line[3]
                        data["permutations"] = line[6]
                        data["mantel_r"] = line[4]
                        data["p_value"] = line[5]
                        data["tail_type"] = line[7]
                        data["cdm"] = "partial_matrix"
                    else:
                        data["entries_num"] = line[2]
                        data["permutations"] = line[5]
                        data["mantel_r"] = line[3]
                        data["p_value"] = line[4]
                        data["tail_type"] = line[6]
                    data_list.append(data)
        try:
            collection = self.db["mantel_test_detail"]
            collection.insert_many(data_list)
        except Exception, e:
            self.bind_object.logger.info("导入mantel检验结果数据出错:%s" % e)
            self.bind_object.set_error("导入mantel检验结果数据出错", code="52801501")
        else:
            self.bind_object.logger.info("导入mantel检验结果数据成功")
