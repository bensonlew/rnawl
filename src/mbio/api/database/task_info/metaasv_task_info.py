# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import datetime
from biocluster.api.database.base import Base, report_check
from bson import SON


class MetaasvTaskInfo(Base):
    def __init__(self, bind_object):
        super(MetaasvTaskInfo, self).__init__(bind_object)
        self._project_type = 'metaasv'

    #@report_check
    def add_task_info(self, db_name=None):
        if db_name:
            self._db_name = db_name
        self.bind_object.logger.info("options: {}".format(self.bind_object.sheet.options))
        json_data = [
            ('task_id', self.bind_object.sheet.id),
            ('member_id', self.bind_object.sheet.member_id),
            ('project_sn', self.bind_object.sheet.project_sn),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('member_type', int(self.bind_object.sheet.member_type)),
            ('cmd_id', int(self.bind_object.sheet.cmd_id)),
            ('is_demo', 0),
            ('demo_id', self.bind_object.sheet.id)
            # ('database_type', (self.bind_object.sheet.options)["database_type"]),
            # ('database', (self.bind_object.sheet.options)["database"])
        ]
        #if "pipeline" in self.bind_object.sheet.options:
        #    json_data.append(('pipeline', (self.bind_object.sheet.options)["pipeline"]))
        self.db['sg_task'].insert_one(SON(json_data))
        self.bind_object.logger.info('任务信息导入sg_task成功。')

    def update_sg_task(self, task_id, database_type, database, pipeline=None):
        """
        更新主表信息
        :param database_type: 更新字段database_type
        :param database:更新字段database
        :param task_id:任务
        :return:
        """
        try:
            if pipeline:
                self.db['sg_task'].update_one({"task_id": task_id}, {"$set": {"database_type": database_type,  "database": database, "pipeline":pipeline}})
            else:
                self.db['sg_task'].update_one({"task_id": task_id}, {"$set": {"database_type": database_type, "database": database, "pipeline": "default"}})
        except:
            self.bind_object.set_error('更新任务失败！')

    def add_sample_numbers(self, task_id, sample_numbers):
        """
        更新主表信息
        :param database_type: 更新字段database_type
        :param database:更新字段database
        :param task_id:任务
        :return:
        """
        try:

            self.db['sg_task'].update_one({"task_id": task_id}, {"$set": {"sample_numbers": sample_numbers}})
        except:
            self.bind_object.set_error('"sample_numbers导入sg_task表格信息出错！')