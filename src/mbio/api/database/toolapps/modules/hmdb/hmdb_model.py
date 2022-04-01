# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
from biocluster.api.database.base import Base, report_check
import datetime
import os
import json
from bson import ObjectId
from bson.son import SON
from biocluster.config import Config

class HmdbModel(Base):
    def __init__(self, bind_object):
        super(HmdbModel, self).__init__(bind_object)
        sanger_type, sanger_path = self.bind_object._sheet.output.split(':')
        sanger_prefix = Config().get_netdata_config(sanger_type)
        self.work_dir = os.path.join(sanger_prefix[sanger_type + '_path'], sanger_path)
        # self.work_dir = self.bind_object.work_dir
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
        self.insert_table(self.work_dir + '/model_result.xls', '预测结果表')
        # self.insert_detail(main_id, self.work_dir + '/annotation_result.xls')

    def insert_table(self, path, name):
        self.bind_object.logger.info('开始导入table表')
        params, risk = self.get_info()
        risk_prob = ""
        prof_list = []
        with open(path,"r") as f1:
            files = f1.readlines()
            for line in files:
                line = line.strip().split("\t")
                if line[0] == "risk_prob":
                    risk_prob = float(line[1])
                else:
                    prof_list.append("%.4f" % float(line[1]))
        risk_str = 'medium'
        if risk_prob < risk[0]:
            risk_str = 'low'
        elif risk_prob > risk[1]:
            risk_str = 'high'

        self.db['hmd_model'].insert_one(SON(
            project_sn=self.bind_object.sheet.project_sn,
            task_id=self.bind_object.sheet.id,
            name=name,
            status='end',
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            position=path,
            prof=prof_list,
            value=risk_prob,
            member_id=self.bind_object.sheet.member_id,
            is_demo=0,
            params = json.dumps(params, sort_keys=True, separators=(',', ':')),
            risk=risk_str,
        ))
        self.bind_object.logger.info('table表导入结束')

    def get_info(self):
        model = self.bind_object.sheet.option('model')
        model_type = self.bind_object.sheet.option('model_type')
        check_coll = self.db['hmd_model_pic']
        result = check_coll.find_one({"disease": model, "model": model_type})
        if not result:
            raise Exception("数据库hmd_model_pic中不存在满足查询条件的数据")
        thresh_list = result['list1']
        params = {
            "model": model,
            "model_type": model_type
        }
        return params,thresh_list

    def check(self):
        """
        检查文件格式是否正确
        """
        pass