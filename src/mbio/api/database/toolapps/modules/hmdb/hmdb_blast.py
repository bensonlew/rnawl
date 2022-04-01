# -*- coding: utf-8 -*-
# __author__ = 'zhengyuan'
from biocluster.api.database.base import Base, report_check
import datetime
import os
from bson import SON
from biocluster.config import Config

class HmdbBlast(Base):
    def __init__(self, bind_object):
        super(HmdbBlast, self).__init__(bind_object)
        sanger_type, sanger_path = self.bind_object._sheet.output.split(':')
        sanger_prefix = Config().get_netdata_config(sanger_type)
        self.work_dir = os.path.join(sanger_prefix[sanger_type + '_path'], sanger_path)
        # self.work_dir = self.bind_object.work_dir
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
        self.insert_table(self.work_dir + '/annotation_result.xls', 'Blast分析结果表')

    def insert_table(self, path, name):
        self.bind_object.logger.info('开始导入table表')
        insert_data = []
        self.db['table'].insert_one(SON(
            project_sn=self.bind_object.sheet.project_sn,
            task_id=self.bind_object.id,
            name=name,
            status='end',
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            position=path
        ))
        self.bind_object.logger.info('table表导入结束')

    def check(self):
        """
        检查文件格式是否正确
        """
        pass