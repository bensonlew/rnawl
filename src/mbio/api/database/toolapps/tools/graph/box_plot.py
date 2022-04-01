# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'
import json
from biocluster.api.database.base import Base, report_check
import re
import os
import datetime
from bson import SON
from biocluster.config import Config
from collections import defaultdict


class BoxPlot(Base):
    def __init__(self, bind_object):
        super(BoxPlot, self).__init__(bind_object)
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
        self.main_id = self.box_plot_in()
        self.table_ids = self.table_in()
        return self.main_id
        pass

    def box_plot_in(self):
        """
        导入box_plot图相关信息
        """
        self.bind_object.logger.info("开始箱线图导表")
        box_file = self.output_dir + '/boxplot.xls'
        box_group_name = []
        with open(box_file) as f:
            box_plot_id = self.db['box_plot'].insert_one(SON(
                project_sn=self.bind_object.sheet.project_sn,
                task_id=self.bind_object.id,
                name='box_plot',
                desc='箱线图的数据',
                status='faild',
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            )).inserted_id
            lines = f.readlines()
            samples = []
            insert_data = []
            for line in lines[1:]:
                insert_filter = []
                line_split = line.strip().split("\t")
                samples.append(line_split[0])
                if len(line_split) == 6:
                    pass
                else:
                    insert_filter = line_split[6].split(";")
                data_list = [float(line_split[1]), float(line_split[2]), float(line_split[3]), float(line_split[4]), float(line_split[5]), insert_filter]
                data = SON(sample_name=line_split[0], box_id=box_plot_id, box_data=data_list)
                insert_data.append(data)
            self.db['box_plot_detail'].insert_many(insert_data)
        if self.bind_object._task.option("group_table").is_set:
            insert_group_data = []
            if self.bind_object._task.option("sed_group") == '':
                all_group_file = os.listdir(self.output_dir)
                for group_file in all_group_file:
                    group_dict = defaultdict(list)
                    if group_file.startswith("group"):
                        group_file = os.path.join(self.output_dir, group_file)
                        with open(group_file)as fr:
                            lines = fr.readlines()
                            group_name = lines[0].strip().split("\t")[1]
                            box_group_name.append(group_name)
                            for line in lines[1:]:
                                sample = line.strip().split("\t")[0]
                                group = line.strip().split("\t")[1]
                                group_dict[group].append(sample)
                            group_data = SON(box_id=box_plot_id, group=group_name, group_data=group_dict)
                            insert_group_data.append(group_data)
            else:
                select_group_file = [self.output_dir + '/first_group.xls', self.output_dir + '/sed_group.xls']
                for group_file in select_group_file:
                    group_dict = defaultdict(list)
                    if group_file.endswith("group.xls"):
                        group_file = os.path.join(self.output_dir, group_file)
                        with open(group_file)as fr:
                            lines = fr.readlines()
                            group_name = lines[0].strip().split("\t")[1]
                            box_group_name.append(group_name)
                            for line in lines[1:]:
                                sample = line.strip().split("\t")[0]
                                group = line.strip().split("\t")[1]
                                group_dict[group].append(sample)
                            group_data = SON(box_id=box_plot_id, group=group_name, group_data=group_dict)
                            print group_data
                            insert_group_data.append(group_data)
            self.db['box_group'].insert_many(insert_group_data)
        self.db['box_plot'].update_one({'_id': box_plot_id}, {'$set': {'status': 'end', 'attrs': samples, 'group_name':box_group_name}})
        return box_plot_id

    def table_in(self):
        """
        导入表格相关信息
        """
        self.bind_object.logger.info("开始表格导表")
        value_table = self.insert_table(self.output_dir + '/boxplot.xls', '箱线图数据表', '箱线图数据表格')
        return [value_table]

    def insert_table(self, fp, name, desc):
        insert_data = []
        attr = ["#name"]
        table_id = self.db['table'].insert_one(SON(
            project_sn=self.bind_object.sheet.project_sn,
            task_id=self.bind_object.id,
            name=name,
            desc=desc,
            status='failed',
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        )).inserted_id
        with open(fp) as f:
            lines = f.readlines()
            columns = lines[0].strip().split("\t")
            for i in columns:
                attr.append(i)
            for line2 in lines[1:]:
                line_split = line2.strip().split('\t')
                data = dict(zip(attr, line_split))
                data['table_id'] = table_id
                insert_data.append(data)
            self.db['table_detail'].insert_many(insert_data)
            self.db['table'].update_one({'_id': table_id}, {'$set': {'status': 'end', 'attrs': columns}})

        return table_id

    def check(self):
        """
        检查文件格式是否正确
        """
        pass
