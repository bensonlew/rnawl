# -*- coding: utf-8 -*-
# __author__ = 'qndanhua'
# import json
from biocluster.api.database.base import Base, report_check
# import re
import os
import datetime
from bson import SON
from biocluster.config import Config


class Chord(Base):
    def __init__(self, bind_object):
        super(Chord, self).__init__(bind_object)
        self.output_dir = self.bind_object.output_dir
        self.work_dir = self.bind_object.work_dir
        if Config().MONGODB == 'sanger':
            self._db_name = 'toolapps'
        else:
            self._db_name = 'ttoolapps'
        self._project_type = 'toolapps'
        self.check()

    @report_check
    def run(self):
        """
        运行函数
        """
        self.main_id = self.add_chord_table(self.output_dir + "", "chord", "弦图", "chord")
        self.table_ids = self.table_in()
        return self.main_id
        pass

    def add_chord_table(self, fp, name, desc, col):

        with open(fp) as f:
            columns = f.readline().strip().split('\t')
            main_data = SON(
                project_sn=self.bind_object.sheet.project_sn,
                task_id=self.bind_object.id,
                name=name,
                attrs=columns,
                desc=desc,
                status='end',
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            )
            table_id = self.db[col].insert_one(main_data).inserted_id

            insert_data = []
            print columns
            for line in f:
                line_split = line.strip().split('\t')
                update_data = dict(zip(columns, line_split[1:]))
                data = {
                    "column_name": line_split[0],
                }
                data.update(update_data)
                data['table_id'] = table_id
                insert_data.append(data)
            self.db[col + '_detail'].insert_many(insert_data)
        return table_id

    def table_in(self):
        """
        导入表格相关信息
        """
        value_table = self.add_chord_table(self.output_dir + '/final_value.xls', '弦图数据', '弦图数据表格', "table")
        return [value_table]
    #
    # def insert_table(self, fp, name, desc):
    #     with open(fp) as f:
    #         lines = f.readlines()
    #         columns = lines[0].strip().split("\t")
    #         insert_data = []
    #         table_id = self.db['table'].insert_one(SON(
    #             project_sn=self.bind_object.sheet.project_sn,
    #             task_id=self.bind_object.id,
    #             name=name,
    #             attrs=columns,
    #             desc=desc,
    #             status='end',
    #             created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    #         )).inserted_id
    #     with open(fp) as f:
    #         lines = f.readlines()
    #         columns = lines[0].strip().split("\t")
    #         for line2 in lines[1:]:
    #             line_split = line2.strip().split('\t')
    #             data = dict(zip(columns, line_split))
    #             data['table_id'] = table_id
    #             insert_data.append(data)
    #         self.db['table_detail'].insert_many(insert_data)
    #     return table_id

    def check(self):
        """
        检查文件格式是否正确
        """
        pass
