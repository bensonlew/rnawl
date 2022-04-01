# -*- coding: utf-8 -*-
from biocluster.api.database.base import Base,report_check
from bson import ObjectId
from bson.son import SON
from mbio.packages.metagenomic.id_convert import name2id
import datetime
import json
import types

class MgAnnoTtss(Base):
    def __init__(self, bind_object):
        super(MgAnnoTtss, self).__init__(bind_object)
        self._project_type = "metagenomic"
        self.sample_2_id = ""  # name2id(self.bind_object.sheet.id, type="task")

    @report_check
    def add_anno_ttss(self, params, anno_file, name="TTSS_Origin"):
        """
        params根据接口设计调整
        geneset_id同其他主表
        anno_file路径同其他主表
        """
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "desc": "TTSS效应蛋白预测",
            "created_ts": created_ts,
            "status": "end",
            "name": name,
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "anno_file": anno_file,
            "is_origin": 1
        }
        collection = self.db["anno_ttss"]
        main_id = collection.insert_one(insert_data).inserted_id
        return main_id

    @report_check
    def add_anno_ttss_stat(self, main_id, table, type="best_hit"):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error("main_id必须为ObjectId对象或其对应的字符串！", code="52803901")
        data_list = list()
        result = self.db["anno_ttss"].find_one({'_id': main_id})
        self.sample_2_id = name2id(result["task_id"], type="task")
        with open(table, "r") as f1:
            file = f1.readlines()
            for line in file[1:]:
                line = line.strip().split("\t")
                data = [
                    ("ttss_id", main_id),
                    ("specimen", self.sample_2_id.get(line[0], line[0])),
                    ("type", type),
                    ("total", float(line[1])),
                    ("is_sec", float(line[2])),
                    ("not_sec", float(line[3])),
                    ("percent", float(line[4]))
                ]
                data = SON(data)
                data_list.append(data)
        detail_coll = self.db["anno_ttss_stat"]
        try:
            detail_coll.insert_many(data_list)
        except Exception,e:
            self.bind_object.set_error( "anno_ttss_stat error: %s" , variables=( e), code="52803902")
        else:
            self.bind_object.logger.info("导入anno_ttss_stat成功")

    @report_check
    def add_anno_ttss_tax(self, main_id, table, type="best_hit"):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error("main_id必须为ObjectId对象或其对应的字符串！", code="52803903")
        data_list = list()
        with open(table, "r") as f1:
            file = f1.readlines()
            for line in file[1:]:
                line = line.strip().split("\t")
                data = [
                    ("ttss_id", main_id),
                    ("type", type),
                    ("tax", line[0]),
                    ("sec_num", float(line[1])),
                    ("total_num", float(line[2])),
                    ("odds", float(line[3])),
                    ("pvalue", float(line[4])),
                    ("correct", float(line[5]))
                ]
                data = SON(data)
                data_list.append(data)
        detail_coll = self.db["anno_ttss_tax"]
        try:
            detail_coll.insert_many(data_list)
            self.db['anno_ttss'].update_one({'_id': main_id}, {'$set':{'main_id':main_id}})
        except Exception,e:
            self.bind_object.set_error("anno_ttss_tax error: %s" , variables=( e), code="52803904")
        else:
            self.bind_object.logger.info("导入anno_ttss_tax成功")

    def export_test(self):
        params = {
            "biandecanshu": "1"
        }
        anno_file = "s3://biandelujing"
        main_id = self.add_anno_ttss(params, anno_file)
        self.add_anno_ttss_stat(main_id, "tmpsummary")
        self.add_anno_ttss_tax(main_id, "ttss_tax", type="best_hit")
        self.add_anno_ttss_tax(main_id, "ttss_tax2", type="lca")