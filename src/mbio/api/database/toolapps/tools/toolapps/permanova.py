# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
import json
from biocluster.api.database.base import Base, report_check
import re
import datetime
from bson import SON
from biocluster.config import Config
import os


class Permanova(Base):
    def __init__(self, bind_object):
        super(Permanova, self).__init__(bind_object)
        self.output_dir = self.bind_object.output_dir ##结果目录
        self.work_dir = self.bind_object.work_dir ##项目目录
        if Config().MONGODB == 'sanger':   ##定义小工具的mongo数据库名称
            self._db_name = 'toolapps'
        else:
            self._db_name = 'ttoolapps'
        self.check()

    @report_check
    def run(self):
        """
        运行函数
        """
        self.table_ids = self.table_in()
        return self.table_ids

    def table_in(self):
        """
		导入表格相关信息
		"""
        permanova = self.insert_table(self.output_dir + '/permanova.xls', 'permanova结果表',
                                     'permanova分析结果数据')
        return permanova

    def insert_table(self, fp, name, desc):
        with open(fp) as f:
            columns = f.readline().strip().split('\t')
            insert_data = []
            self.bind_object.logger.info("permanova主表的permanova写入")
            table_id = self.db['table'].insert_one(SON(
                project_sn=self.bind_object.sheet.project_sn,
                task_id=self.bind_object.id,
                name=name,
                attrs=columns,
                desc=desc,
                status='fail',
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )).inserted_id
            for line in f:
                line = line.strip().split("\t")
                data = [("table_id", table_id), ("Characteristics", line[0]), ("SumsOfSqs", line[1]), ("MeanSqs", line[2]),
                        ("F_Model", line[3]), ("R2", line[4]), ("p_value", line[5]),("p_adjust",  line[6])]
                data_son = SON(data)
                insert_data.append(data_son)
        self.db['table_detail'].insert_many(insert_data)
        try:
            collection = self.db["table"]
            collection.update_one({"_id": table_id}, {"$set": {"status":'end'}})
        except Exception, e:
            self.bind_object.logger.error("导入permanova信息出错:%s" % e)
        else:
            self.bind_object.logger.info("导入permanova信息成功!")
        return table_id


    def check(self):
        """
        检查文件格式是否正确
        """
        pass
