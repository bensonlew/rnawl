# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
import json
from biocluster.api.database.base import Base, report_check
import re
import datetime
from bson import SON
from biocluster.config import Config


class Circos(Base):
    def __init__(self, bind_object):
        super(Circos, self).__init__(bind_object)
        self.output_dir = self.bind_object.output_dir
        self.work_dir = self.bind_object.work_dir
        if Config().MONGODB == 'sanger':
            self._db_name = 'toolapps'
        else:
            self._db_name = 'ttoolapps'
        self._project_type = 'toolapps'
        # self.check()

    @report_check
    def run(self):
        """运行函数"""
        if self.bind_object._task.option("group_table").is_set:
            group_name = self.bind_object._task.option("group_table").prop["group_scheme"]
            self.bind_object.logger.info(group_name)
            for g in group_name:
                fp = self.output_dir + "/{}_result_data".format(g)
                self.main_id = self.circos_in(fp, g)
                self.table_ids = self.table_in(fp, g)
        else:
            fp = self.output_dir + "/all_result_data"
            self.main_id = self.circos_in(fp)
            self.table_ids = self.table_in(fp)

    def table_in(self, fp, group=None):
        """导入表格信息"""
        self.insert_table(fp, 'circos图结果表', '画circos图时使用的数据', group_name=group)

    def insert_table(self, fp, name, desc, group_name=None):
        with open(fp) as f:
            columns = f.readline().strip().split('\t')
            insert_data = []
            table_id = self.db['table'].insert_one(SON(
                project_sn=self.bind_object.sheet.project_sn,
                task_id=self.bind_object.id,
                name=name,
                attrs=columns,
                desc=desc,
                status='end',
                group_name=group_name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )).inserted_id
            for line in f:
                line_split = line.strip().split('\t')
                data = dict(zip(columns, line_split))
                data['table_id'] = table_id
                insert_data.append(data)
            self.db['table_detail'].insert_many(insert_data)
        return table_id

    def circos_in(self, fp, group_name=None):
        """导入circos画图数据"""
        with open(fp) as f:
            self.bind_object.logger.info("circos主表导入")
            head = f.next().strip('\r\n')  # windows换行符
            head = re.split('\t', head)
            new_head = head[1:]
            circos_id = self.db['circos'].insert_one(SON(
                project_sn=self.bind_object.sheet.project_sn,
                task_id=self.bind_object.id,
                name='circos',
                desc='弦图',
                status='failed',
                attrs=new_head,
                group_name=group_name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            )).inserted_id
            insert_data = []
            for line in f:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                sample_num = line[1:]
                otu_detail = dict()
                otu_detail['circos_id'] = circos_id
                otu_detail['name'] = line[0]
                for i in range(0, len(sample_num)):
                    otu_detail[new_head[i]] = float(sample_num[i])  # 保证画图时取到的数据是数值型
                insert_data.append(otu_detail)
            self.bind_object.logger.info("circos画图数据导入开始")
            try:
                self.db['circos_detail'].insert_many(insert_data)
                self.db['circos'].update_one({'_id': circos_id}, {'$set': {'status': 'end'}})
            except Exception as e:
                self.bind_object.logger.info("circos画图数据导入出错{}".format(e))
                raise Exception("circos画图数据导入出错{}".format(e))
            else:
                self.bind_object.logger.info("circos画图数据导入完成")
