# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
# last modified by guhaidong (add member_type cmd_id)
import datetime
from biocluster.api.database.base import Base, report_check
# from biocluster.config import Config
from bson import SON


class MgTaskInfo(Base):
    def __init__(self, bind_object):
        super(MgTaskInfo, self).__init__(bind_object)
        self._project_type = 'metagenomic'
        #self._db_name = Config().MONGODB + "_metagenomic"

    #@report_check
    def add_task_info(self, db_name=None):
        if db_name:
            self._db_name = db_name
        json_data = [
            ('task_id', self.bind_object.sheet.id),
            ('member_id', self.bind_object.sheet.member_id),
            ('project_sn', self.bind_object.sheet.project_sn),
            ('created_ts', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            ('member_type', int(self.bind_object.sheet.member_type)),
            ('cmd_id', int(self.bind_object.sheet.cmd_id)),
            ('is_demo', 0),
            ('demo_id', self.bind_object.sheet.id)
        ]
        self.db['sg_task'].insert_one(SON(json_data))
        self.bind_object.logger.info('任务信息导入sg_task成功。')
