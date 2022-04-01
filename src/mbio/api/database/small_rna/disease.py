# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'
# last_modify: liubinxu 20181121
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
from biocluster.config import Config
import pandas as pd
import datetime
import json
import unittest
import os

from mbio.api.database.small_rna.api_base import ApiBase


class Disease(ApiBase):
    def __init__(self, bind_object):
        super(Disease, self).__init__(bind_object)

    def import_disease(self, species=None, mir2disease=None, hmdd=None):
        task_id = self.bind_object.sheet.id
        if species != "Homo_sapiens":
            self.bind_object.logger.info("物种不是Homo_sapiens")
            return

        self.remove_table_by_main_record(main_table='sg_disease', task_id=task_id, detail_table=['sg_disease_mir2', 'sg_disease_hmdd'], detail_table_key='disease_id')
        main_id = self.add_disease_main()
        self.add_minet_mir2(main_id, mir2disease)
        self.add_minet_hmdd(main_id, hmdd)
        self.db['sg_disease'].update({"_id": main_id},
                                   {"$set": {"status": "end", "main_id": main_id}})

    @report_check
    def add_disease_main(self):
        data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "status": "start",
            "name": 'miRNA_disease_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3],
            "desc": "miRNA disease",
            "created_ts": datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        }
        col = self.db["sg_disease"]
        main_id = col.insert_one(data).inserted_id
        self.bind_object.logger.info("主表创建成功")
        return main_id

    @report_check
    def add_minet_mir2(self, main_id,  mir2disease):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, StringTypes):
                main_id = ObjectId(main_id)
            else:
                raise Exception("minet_id必须为ObjectId对象或者其对应的字符串！")
        data_list = []
        with open(mir2disease, "rb") as mir2disease_f:
            # head_line = mir2disease_f.readline()
            #  header = head_line.strip("\n").split("\t")
            header = [
                "mirna",
                "disease",
                "style",
                "method",
                "date",
                "description",
            ]
            for line in mir2disease_f:
                cols = line.strip("\n").split("\t")
                data = zip(header, cols)
                data.append(('disease_id', main_id))
                data = SON(data)
                data_list.append(data)

        try:
            self.create_db_table('sg_disease_mir2', data_list)
        except Exception as e:
            self.bind_object.set_error("导入MIR2DISEASE信息出错!")
        else:
            self.bind_object.logger.info("导入MIR2DISEASE信息成功!")

    @report_check
    def add_minet_hmdd(self, main_id, hmdd):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, StringTypes):
                main_id = ObjectId(main_id)
            else:
                raise Exception("main_id必须为ObjectId对象或者其对应的字符串！")
        data_list = []
        with open(hmdd, "rb") as hmdd_f:
            head_line = hmdd_f.readline()
            header = head_line.strip("\n").split("\t")
            header[1] = "mirna"
            for line in hmdd_f:
                cols = line.strip("\n").split("\t")
                data = zip(header, cols)
                data.append(('disease_id', main_id))
                data = SON(data)
                data_list.append(data)

        try:
            self.create_db_table('sg_disease_hmdd', data_list)
        except Exception as e:
            self.bind_object.set_error("导入HMDD出错!")
        else:
            self.bind_object.logger.info("导入HMDD信息成功!")


if __name__ == '__main__':

    class TestFunction(unittest.TestCase):
        """
        测试导表函数
        """
        def test_mongo(test):
            from mbio.workflows.small_rna.small_rna_test_api import SmallRnaTestApiWorkflow
            from biocluster.wsheet import Sheet
            import random
            data = {
                "id": "small_rna",
                #+ str(random.randint(1,10000)),
                #"id": "denovo_rna_v2",
                "project_sn": "small_rna",
                #+ str(random.randint(1,10000)),
                "type": "workflow",
                "name": "small_rna.small_rna_test_api",
                "options": {

                }
            }
            wsheet = Sheet(data=data)
            wf = SmallRnaTestApiWorkflow(wsheet)

            wf.IMPORT_REPORT_DATA = True
            wf.IMPORT_REPORT_AFTER_END = False
            wf.test_api = wf.api.api("small_rna.disease")

            wf.test_api.import_disease("Homo_sapiens", "/mnt/lustre/users/sanger/app/database/miR2Disease/AllEntries.txt" ,"/mnt/lustre/users/sanger/app/database/miR2Disease/alldata.txt")
    unittest.main()
