# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'
import json
from biocluster.api.database.base import Base, report_check
import re
import os
import datetime
from bson import SON
from biocluster.config import Config


class SimplePie(Base):
    def __init__(self, bind_object):
        super(SimplePie, self).__init__(bind_object)
        self.output_dir = self.bind_object.output_dir
        self.work_dir = self.bind_object._task.work_dir
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
        if self.bind_object._task.option("group_table").is_set and self.bind_object._task.option("calculation") != 'none':
            for i in self.bind_object._task.option("group_table").prop['group_scheme']:
                self.main_id = self.simple_pie_in(i)
                self.table_ids = self.table_in(i)
        else:
            self.main_id = self.simple_pie_in()
            self.table_ids = self.table_in()
        return self.main_id
        pass

    def simple_pie_in(self, group=None):
        """
        导入simple_pie图相关信息
        """
        self.bind_object.logger.info("开始饼图导表")
        if self.bind_object._task.option("group_table").is_set and self.bind_object._task.option("calculation") != 'none':
            pie = self.insert_pie(self.output_dir + '/' + group + '_matrix_pie.xls', 'pie', '一个样本多个值的饼图表格', group)

        else:
            pie = self.insert_pie(self.output_dir + '/matrix_pie.xls', 'pie', '一个样本多个值的饼图表格')
        return [pie]

    def insert_pie(self, fp, name, desc, group_name=None):
        with open(fp) as f:
            simple_pie_id = self.db['pie'].insert_one(SON(
                project_sn=self.bind_object.sheet.project_sn,
                task_id=self.bind_object.id,
                name=name,
                desc=desc,
                status='faild',
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                group_name=group_name,
            )).inserted_id
            lines = f.readlines()
            lines = [line for line in lines if (line != "\r\n") and (line != "\n")]
            lines = [line for line in lines if not re.search(r"^(\s*\t+?)\s*\t*\n*", line)]
            first_line = lines[0].strip().split("\t")
            xAxis = []
            for i in first_line[1:]:
                xAxis.append(i)
            samples = []
            insert_data = []
            for line in lines[1:]:
                sample_data = []
                line_split = line.strip().split("\t")
                samples.append(line_split[0])
                for i in range(len(first_line)-1):
                        dic = dict()
                        dic["name"] = first_line[i+1]
                        dic["value"] = float(line_split[i+1])
                        sample_data.append(dic)
                data = SON(sample_name=line_split[0], pie_id=simple_pie_id, value=sample_data)
                insert_data.append(data)
            self.db['pie_detail'].insert_many(insert_data)
            self.db['pie'].update_one({'_id': simple_pie_id},
                                        {'$set': {'status': 'end', 'attrs': samples, 'categories': xAxis}})
            return simple_pie_id

    def table_in(self, group=None):
        """
        导入表格相关信息
        """
        self.bind_object.logger.info("开始表格导表")
        if self.bind_object._task.option("group_table").is_set and self.bind_object._task.option("calculation") != 'none':
            value_table = self.insert_table(self.output_dir + '/' + group + '_final_value.xls', '饼图数据表', '饼图的数据表格',
                                            group)
        else:
            value_table = self.insert_table(self.output_dir + '/final_value.xls', '饼图数据表', '饼图的数据表格')
        return [value_table]

    def insert_table(self, fp, name, desc, group_name=None):
        with open(fp) as f:
            lines = f.readlines()
            columns = lines[0].strip().split("\t")
            insert_data = []
            table_id = self.db['table'].insert_one(SON(
                project_sn=self.bind_object.sheet.project_sn,
                task_id=self.bind_object.id,
                name=name,
                attrs=columns,
                desc=desc,
                status='end',
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                group_name=group_name,
            )).inserted_id
        with open(fp) as f:
            lines = f.readlines()
            columns = lines[0].strip().split("\t")
            for line2 in lines[1:]:
                line_split = line2.strip().split('\t')
                data = dict(zip(columns, line_split))
                data['table_id'] = table_id
                insert_data.append(data)
            self.db['table_detail'].insert_many(insert_data)
        return table_id

    def check(self):
        """
        检查文件格式是否正确
        """
        pass
