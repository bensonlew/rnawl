# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
#20191104

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re
from types import StringTypes
from biocluster.config import Config
from bson.son import SON


class PanVenn(Base):
    def __init__(self, bind_object):
        super(PanVenn, self).__init__(bind_object)
        self._project_type = "bac_comparative"
        # self.id = 'tsg_123'
        # self.project_sn = '188_5b5acb3018'

    @report_check
    def add_venn(self, params=None, name=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        # task_id = self.id
        # project_sn = self.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "Venn图主表",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "created_ts": created_ts,
            "name": name if name else "Pan_Venn_Origin",
        }
        collection = self.db["pan_venn"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': inserted_id},{'$set': {'main_id':inserted_id}})
        return inserted_id

    def add_venn_detail(self, inserted_id, venn_path):
        if not isinstance(inserted_id, ObjectId):
            new_inserted_id = ObjectId(inserted_id)
        else:
            new_inserted_id = inserted_id
        data_list = []
        with open(venn_path, "r") as f:
            lines = f.readlines()
            for line in lines:
                if (re.search(r"Group_name", line)) or (re.search(r"Sample_name", line)):
                    pass
                else:
                    line = line.strip().split("\t")
                    index_name = line[0].split(" only")[0].strip()
                    gene_num = line[1]
                    try:
                        gene_list = line[2]
                    except:
                        gene_list = "-"
                    data = {
                        "pan_id": new_inserted_id,
                        "label": index_name,
                        "cluster_list": gene_list,
                        "gene_num": gene_num,
                        }
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            collection = self.db["pan_venn_stat"]
            collection.insert_many(data_list)
            main_collection = self.db["pan_venn"]
            main_collection.update_one({'_id': new_inserted_id},{'$set': {'main_id':new_inserted_id}})
        except Exception, e:
            self.bind_object.logger.info("导入%s结果表出错:%s" % (venn_path, e))