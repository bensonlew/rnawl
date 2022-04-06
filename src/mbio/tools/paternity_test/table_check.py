# !/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/python
# -*- coding: utf-8 -*-

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import os
import re
import shutil
import MySQLdb


class TableCheckAgent(Agent):
    """
    该tool用于将连接limus系统，下载客户信息表，并同步到monogo中，同时合并上机表
    version v1.0
    author: hongdongxuan
    last_modify: 20171013
    """
    def __init__(self, parent):
        super(TableCheckAgent, self).__init__(parent)
        options = [
            {"name": "table_up", "type": "infile", 'format': 'paternity_test.split_check'},  # 输入线上上机表
            {"name": "table_down", "type": "infile", "format": "paternity_test.split_check"}  # 输入线下上机表
        ]
        self.add_option(options)
        self.step.add_steps("table_check")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.table_check.start()
        self.step.update()

    def stepfinish(self):
        self.step.table_check.finish()
        self.step.update()


    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("table_up").is_set:
            raise OptionError("必须输入线上上机表")
        if not self.option("table_down").is_set:
            raise OptionError("必须输入线下上机表")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(TableCheckAgent, self).end()


class TableCheckTool(Tool):
    """
    运行获取客户信息，然后导入到mongo表中
    """
    def __init__(self, config):
        super(TableCheckTool, self).__init__(config)
        self._version = '1.0.1'
        self.DB_HOST = "172.16.101.202"
        self.DB_USER = 'mjlimsbak'
        self.DB_PASSWD = 'Q6k9pU9ZN5FTow'
        self.DB_NAME = 'mjlims'

    def update_pt_customer(self):
        """
        用于更新亲子鉴定的客户信息表，到sg_guest_info
                case_name : data[0]  # 家系编号  WQ17103749
                create_time = data[1]  # 创建时间
                ask_person = data[2]  # 申请人 张文
                sample_name = data[3]  # 母本与父本名称 尚雪
                sample_type = data[4]  # 母本与父本样本类型，全血,白细胞等
                sample_number = data[5]  # 父本与母本编号 BB1710060001
                sample_id = data[6]  # 母本与父本id M-1
                ask_time = data[7]  # 申请时间
                accept_time = data[8]  # 受理时间
                isurgent = data[9]  # 是否加急
                gestation_week = data[10]  # 孕周
                company_name = data[11]  # 代理公司
                contacts_people = data[12]  # 代理公司联系人
        :return:
        """
        self.option("table_up").get_info()
        wq_sample_list = self.option("table_up").prop['wq_sample_list']
        self.logger.info("亲子数据列表{}".format(wq_sample_list))
        self._db_client = MySQLdb.connect(host=self.DB_HOST, user=self.DB_USER, passwd=self.DB_PASSWD, db=self.DB_NAME,
                                          charset='utf8')
        sql = "SELECT i.note, o.create_date, o.name, i.patient_name, st.name, i.code, i.parented, o.sjrq, " \
              "i.accept_date, o.isurgent, o.gestational_weeks, t.name, t.principal FROM sample_info i " \
              "LEFT JOIN sample_order o ON i.sample_order = o.id LEFT JOIN primary_task t ON t.id = o.advance " \
              "LEFT JOIN dic_sample_type st ON st.id = i.sample_kind WHERE o.id IS NOT NULL AND i.note LIKE 'WQ%'" \
              " and i.note in {} ORDER BY o.create_date DESC".format(tuple(wq_sample_list))
        cursor = self._db_client.cursor()
        cursor.execute(sql)
        results = cursor.fetchall()
        data_list_insert = []
        data_list_update = []
        for data in results:
            if self.api.paternity_test_v2.find_sample_id(data[0] + data[6]):
                update_data = {
                    "case_name": data[0],
                    "create_time": data[1],
                    "ask_person": data[2],
                    "sample_name": data[3],
                    "sample_type": data[4],
                    "sample_number": data[5],
                    "sample_id": data[6],
                    "ask_time": data[7],
                    "accept_time": data[8],
                    "isurgent": data[9],
                    "gestation_week": data[10],
                    "company_name": data[11],
                    "contacts_people": data[12],
                    'insert_time': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                }
                data_list_update.append(update_data)
            else:
                insert_data = {
                    "case_name": data[0],
                    "create_time": data[1],
                    "ask_person": data[2],
                    "sample_name": data[3],
                    "sample_type": data[4],
                    "sample_number": data[5],
                    "sample_id": data[6],
                    "ask_time": data[7],
                    "accept_time": data[8],
                    "isurgent": data[9],
                    "gestation_week": data[10],
                    "company_name": data[11],
                    "contacts_people": data[12],
                    "name": data[0] + data[6],
                    'insert_time': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                }
                data_list_insert.appned(insert_data)
        self._db_client.close()
        self.logger.info(data_list_insert)
        self.logger.info("- - - - - - - - -")
        self.logger.info(data_list_update)
        self.api.paternity_test_v2.add_pt_guest(data_list_insert, "insert")
        self.api.paternity_test_v2.add_pt_guest(data_list_update, "update")

    def update_nipt_customer(self):
        """
        用于更新产筛的客户信息表，到sg_guest_info
        :return:
        """
        pass

    def run(self):
        super(TableCheckTool, self).run()
        self.update_pt_customer()
        self.end()
