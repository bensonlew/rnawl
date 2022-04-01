# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
from biocluster.api.database.base import Base, report_check
import datetime
import os
import json
import types
from bson import ObjectId
from bson.son import SON
from biocluster.config import Config
import pandas as pd

class Train(Base):
    def __init__(self, bind_object):
        super(Train, self).__init__(bind_object)
        # sanger_type, sanger_path = self.bind_object._sheet.output.split(':')
        # sanger_prefix = Config().get_netdata_config(sanger_type)
        # self.work_dir = os.path.join(sanger_prefix[sanger_type + '_path'], sanger_path)
        self.work_dir = self.bind_object.work_dir
        if Config().MONGODB == 'sanger':
            self._db_name = 'hmdb'
        else:
            self._db_name = 'hmdb'
        self._project_type = 'hmdb'
        self.check()

    @report_check
    def run(self):
        """
        运行函数
        """
        disease = self.bind_object.sheet.option("disease")
        model_type = self.bind_object.sheet.option("model_type")
        # table_path = os.path.join(self.work_dir, disease + "_" + model_type + "_model")
        table_path = os.path.join(self.work_dir, "Train", disease + "_" + model_type + "_model")
        # table_path = "" # self.bind_object._tools_report_data[0].option("out").path
        self.bind_object.logger.info(table_path)
        main_id = self.insert_table(table_path, '训练模型表')
        self.insert_pic(main_id, table_path + ".draw", "mt_pic")
        self.insert_pic(main_id, table_path + ".auc", "mt_auc")

    def insert_table(self, path, name):
        self.bind_object.logger.info('开始导入table表')
        train_file = self.bind_object.sheet.option("train_file")
        disease = self.bind_object.sheet.option("disease")
        model_type = self.bind_object.sheet.option("model_type")
        params = {
            "train_file": train_file,
            "disease": disease,
            "model_type": model_type
        }
        class_0, class_1, class_total, auc = self.get_report(path + ".report")
        main_id = self.db['mt'].insert_one(SON(
            project_sn=self.bind_object.sheet.project_sn,
            task_id=self.bind_object.sheet.id,
            name=name,
            status='end',
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            position=os.path.join(self.bind_object.sheet.output, os.path.basename(path)),
            member_id=self.bind_object.sheet.member_id,
            is_demo=0,
            params=json.dumps(params, sort_keys=True, separators=(',', ':')),
            class_0=class_0,
            class_1=class_1,
            class_total=class_total,
            auc=auc
        )).inserted_id
        self.bind_object.logger.info('table表导入结束')
        return main_id

    def get_report(self, path):
        with open(path, "r") as file:
            file.readline()
            file.readline()
            class_0 = file.readline().strip().split()[-4:-1]
            class_1 = file.readline().strip().split()[-4:-1]
            file.readline()
            total = file.readline().strip().split()[-4:-1]
            auc = file.readline().strip().split()[-1]
        return class_0, class_1, total, auc


    def insert_pic(self, main_id, path, db_name):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                raise Exception('main_id必须为ObjectId对象或其对应的字符串！')
        data = pd.read_table(path)
        data["model_id"] = main_id
        insert_data = data.to_dict(orient="records")
        detail_coll = self.db[db_name]
        try:
            detail_coll.insert_many(insert_data)
        except Exception,e:
            self.bind_object.logger.info("%s 表导入失败" % db_name)
            self.bind_object.logger.info(e)
        else:
            self.bind_object.logger.info("%s 表导入成功" % db_name)

    def check(self):
        """
        检查文件格式是否正确
        """
        pass