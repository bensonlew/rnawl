# -*- coding: utf-8 -*-
# __author__ = 'moli.zhou'
import datetime
import random

from biocluster.api.database.base import Base, report_check
import re
from biocluster.config import Config
from pymongo import MongoClient
import gridfs
from mainapp.libs.param_pack import param_pack
from bson import ObjectId
import xlrd

class NiptFlow(Base):
    '''
    将前端需要调用的结果文件导入mongo数据库，之结果保存的tsanger collection
    '''
    def __init__(self, bind_object):
        super(NiptFlow, self).__init__(bind_object)
        self.mongo_client = Config().mongo_client
        self.database = self.mongo_client[Config().MONGODB+'_nipt']

    def nipt_customer(self,file):
        bk = xlrd.open_workbook(file)
        sh = bk.sheet_by_name(u'\u65e0\u521b\u4ea7\u7b5b')
        nrows = sh.nrows

        insert = list()
        # 获取各行数据
        for i in range(4, nrows):
            row_data = sh.row_values(i)
            if i == 4:
                report_num_index = row_data.index(
                    u'\u9879\u76ee\u7f16\u53f7\uff08\u62a5\u544a\u7f16\u53f7\uff09')  # 报告编号
                sample_date_index = row_data.index(u'\u91c7\u6837\u65e5\u671f')  # 采样日期
                patient_name_index = row_data.index(u'\u60a3\u8005\u59d3\u540d\u6837\u672c\u540d\u79f0')  # 患者姓名
                accpeted_date_index = row_data.index(u'\u6536\u6837\u65e5\u671f')  # 收样日期
                number_index = row_data.index(u'\u7f8e\u5409\u6761\u5f62\u7801')  # 样本编号
                register_number_index = row_data.index(u'\u4f4f\u9662/\u95e8\u8bca\u53f7')  # 住院号
                gestation_index = row_data.index(u'\u5b55\u4ea7\u53f2')  # 孕产史
                final_period_index = row_data.index(u'\u672b\u6b21\u6708\u7ecf')  # 末次月经
                gestation_week_index = row_data.index(u'\u5b55\u5468')  # 孕周
                pregnancy_index = row_data.index(u'\u5355\u80ce/\u53cc\u80ce')  # 单双胎
                IVFET_index = row_data.index(u'IVF-ET\u598a\u5a20')  # IVF-ET妊娠
                hospital_index = row_data.index(u'\u9001\u68c0\u5355\u4f4d/\u533b\u9662')  # 送检单位、医院
                doctor_index = row_data.index(u'\u9001\u68c0\u533b\u751f')  # 送检医生
                tel_index = row_data.index(u'\u60a3\u8005\u8054\u7cfb\u7535\u8bdd')  # 患者联系方式
                status_index = row_data.index(u'\u6807\u672c\u72b6\u6001\u5f02\u5e38')  # 标本状态异常
                age_index = row_data.index(u'\u5e74\u9f84')  # 年龄
                type_index = row_data.index(u'\u6807\u672c\u7c7b\u578b')  # 样本类型
            else:
                print i
                para_list = []
                report_num = row_data[report_num_index]
                para_list.append(report_num)
                sample_date = row_data[sample_date_index]
                if type(sample_date) == 'float':
                    sample_date_tuple = xlrd.xldate_as_tuple(sample_date, 0)
                    para_list.append(
                        str(sample_date_tuple[0]) + '/' + str(sample_date_tuple[1]) + '/' + str(sample_date_tuple[2]))
                else:
                    para_list.append(sample_date)

                patient_name = row_data[patient_name_index]
                para_list.append(patient_name)
                accpeted_date = row_data[accpeted_date_index]
                
                if type(accpeted_date) == 'float':
                    accpeted_date_tuple = xlrd.xldate_as_tuple(accpeted_date, 0)
                    para_list.append(
                    str(accpeted_date_tuple[0]) + '/' + str(accpeted_date_tuple[1]) + '/' + str(accpeted_date_tuple[2]))
                else:
                    para_list.append(accpeted_date)

                number = row_data[number_index]
                para_list.append(number)
                register_number = row_data[register_number_index]
                para_list.append(register_number)
                gestation = row_data[gestation_index]
                para_list.append(gestation)
                pregnancy = row_data[pregnancy_index]
                para_list.append(pregnancy)
                IVFET = row_data[IVFET_index]
                para_list.append(IVFET)
                hospital = row_data[hospital_index]
                para_list.append(hospital)
                doctor = row_data[doctor_index]
                para_list.append(doctor)
                tel = row_data[tel_index]
                para_list.append(str(tel))
                status = row_data[status_index]
                para_list.append(status)
                final_period = row_data[final_period_index]
                para_list.append(final_period)
                gestation_week = row_data[gestation_week_index]
                para_list.append(gestation_week)

                age = row_data[age_index]
                para_list.append(str(age))
                sample_type = row_data[type_index]
                para_list.append(sample_type)

                collection = self.database['nipt_customer']
                if collection.find_one({"report_num":para_list[0]}):
                    continue
                else:
                    insert_data = {
                        "report_num":para_list[0],
                        "sample_date": para_list[1],
                        "patient_name": patient_name,
                        "accpeted_date": para_list[3],
                        "number": para_list[4],
                        "register_number": para_list[5],
                        "gestation": para_list[6],
                        "pregnancy": para_list[7],
                        "IVFET": para_list[8],
                        "hospital": para_list[9],
                        "doctor": para_list[10],
                        "tel": para_list[11],
                        "status": para_list[12],
                        "final_period": para_list[13],
                        "gestation_week": para_list[14],
                        "age":para_list[15],
                        "sample_type":para_list[16],
                    }
                    insert.append(insert_data)
        try:
            collection = self.database['nipt_customer']
            collection.insert_many(insert)
        except Exception as e:
            raise Exception('插入客户信息表出错：{}'.format(e))
        else:
            self.bind_object.logger.info("插入客户信息表成功")

