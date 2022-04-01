# -*- coding: utf-8 -*-
from biocluster.api.database.base import Base,report_check
from bson import ObjectId
from bson.son import SON
from mbio.packages.metagenomic.id_convert import name2id
import datetime
import json
import types

class MgAnnoSec(Base):
    def __init__(self, bind_object):
        super(MgAnnoSec, self).__init__(bind_object)
        self._project_type = "metagenomic"
        self.sample_2_id = ""  # name2id(self.bind_object.sheet.id, type="task")

    @report_check
    def add_anno_sec(self, params, anno_file, name="AnnoSec_Origin"):
        """
        params根据接口设计调整
        anno_file路径同其他主表
        """
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "desc": "分泌蛋白预测",
            "created_ts": created_ts,
            "status": "end",
            "name": name,
            # "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            "anno_file": anno_file,
            "is_origin": 1
        }
        collection = self.db["anno_sec"]
        main_id = collection.insert_one(insert_data).inserted_id
        params['anno_id'] = str(main_id)
        collection.update_one({"_id": main_id}, {'$set': {'params': json.dumps(params, sort_keys=True, separators=(',', ':'))}})
        return main_id

    @report_check
    def add_anno_sec_stat(self, main_id, table, model):
        """
        分泌蛋白个数统计表
        :param main_id:
        :param table:
        :param model:  物种类型bac_fun bac fun bac_pos bac_neg
        :return:
        """
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error("main_id必须为ObjectId对象或其对应的字符串！", code="52804401")
        data_list = list()
        result = self.db["anno_sec"].find_one({'_id': main_id})
        self.sample_2_id = name2id(result["task_id"], type="task")
        with open(table, "r") as f1:
            file = f1.readlines()
            for line in file[1:]:
                line = line.strip().split("\t")
                data = [
                    ("sec_id", main_id),
                    ("specimen", self.sample_2_id.get(line[0], line[0])),
                    ("total", float(line[1])),
                    ("true", float(line[2])),
                    ("false", float(line[3])),
                    ("percent", float(line[4])),
                    ("model", model)
                ]
                data = SON(data)
                data_list.append(data)
        detail_coll = self.db["anno_sec_stat"]
        try:
            detail_coll.insert_many(data_list)
        except Exception,e:
            self.bind_object.set_error("anno_sec_stat error: %s" , variables=( e), code="52804402")
        else:
            self.bind_object.logger.info("导入anno_sec_stat成功")

    @report_check
    def add_anno_sec_tax(self, main_id, table, model):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error("main_id必须为ObjectId对象或其对应的字符串！", code="52804403")
        data_list = list()
        with open(table, "r") as f1:
            file = f1.readlines()
            for line in file[1:]:
                line = line.strip().split("\t")
                data = [
                    ("sec_id", main_id),
                    ("model", model),
                    ("tax", line[0]),
                    ("sec_num", float(line[1])),
                    ("total_num", float(line[2])),
                    ("pvalue", float(line[4])),
                    ("correct", float(line[5])),
                    ("odds", float(line[3]))
                ]
                data = SON(data)
                data_list.append(data)
        detail_coll = self.db["anno_sec_tax"]
        try:
            detail_coll.insert_many(data_list)
            self.db['anno_sec'].update_one({'_id': main_id}, {'$set':{'main_id':main_id}})
        except Exception,e:
            self.bind_object.set_error("anno_sec_tax error: %s" , variables=( e), code="52804404")
        else:
            self.bind_object.logger.info("导入anno_sec_tax成功")

    @report_check
    def update_anno_file(self, main_id, anno_file):
        main_coll = self.db["anno_sec"]
        anno_real_file = os.path.join(self.bind_object.output, anno_file)
        main_coll.update_one({"_id": main_id}, {'$set': {'anno_file': anno_real_file}})

    def export_test(self):
        params = {
            "biandecanshu": "1"
        }
        anno_file = "s3://biandelujing"
        main_id = self.add_anno_sec(params, anno_file)
        self.add_anno_sec_stat(main_id, "tmpsummary")
        self.add_anno_sec_tax(main_id, "ttss_tax")