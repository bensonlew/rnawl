# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'xueqinwen'

from collections import OrderedDict
from bson import SON
from bson.objectid import ObjectId
from biocluster.config import Config
from types import StringTypes
from io import StringIO
import types
import re
import datetime
import os
import json
import time
# from api_base import ApiBase
from biocluster.api.database.base import Base, report_check

class BacIdentity(Base):
    def __init__(self, bind_object):
        super(BacIdentity, self).__init__(bind_object)
        self._project_type = 'datasplit'

    def update_bac_db_record(self, main_id, s3_upload_dir,final_dir):
        """
        更新主表字段
        """
        update_dict = {
            "s3": s3_upload_dir,
            "final_name":final_dir
        }
        self.db['bac_identity'].update({"_id":self.check_objectid(main_id)},{'$set':update_dict}, upsert =True, multi= True)


    def add_qc_report(self, main_id, output_file):
        """
        导入sangerReport表
        """
        insert_datas = []
        insert_normals = []
        with open(output_file,"r") as out:
            out.readline()
            while 1:
                line = out.readline()
                if not line:
                    break 
                temp = line.rstrip().split("\t")
                insert_data = {
                    "bac_identity_id":self.check_objectid(main_id),
                    "Majorbio-No":temp[0],
                    "Primer":temp[1],
                    "rawSeqLength":temp[2],
                    "trimSeqLength":temp[3],
                    "QC-Identity":temp[4],
                    "Quality-Rank":temp[5],
                    "Date":temp[6]
                }
                insert_normal={
                    "productionstatus":"完成",
                    "productionID":temp[6],
                    "majorbio-no":temp[0],
                    "seqstatus": "测序完成，测序完成",
                    "seqfeature":temp[1] + "," + temp[5],
                    "tag": "OK"
                }
                insert_normals.append(insert_normal)
                insert_datas.append(insert_data)
        self.col_insert_data("bac_normal_sequencing", insert_normals)
        self.col_insert_data("bac_identity_qc_report_detail", insert_datas)


    def add_blastnt_tax(self, main_id, output_file):
        """
        导入blast结果表
        """
        insert_datas = []
        with open(output_file,"r") as out:
            out.readline()
            while 1:
                line = out.readline()
                if not line:
                    break 
                if "NO-TaxID" in line:
                    continue
                temp = line.rstrip().split("\t")
                insert_data = {
                    "bac_identity_id":self.check_objectid(main_id),
                    "Query":temp[0],
                    "Taxonomy":temp[1],
                    "Subject":temp[2],
                    "Identity":temp[3],
                    "Align_length":temp[4],
                    "Mismatches":temp[5],
                    "Gap_open":temp[6],
                    "Q_start":temp[7],
                    "Q_end":temp[8],
                    "S_start":temp[9],
                    "S_end":temp[10],
                    "Evalue":temp[11],
                    "Bit_score":temp[12],
                }
                insert_datas.append(insert_data)
        self.col_insert_data("bac_blastnt_tax_report_detail", insert_datas)

    def col_insert_data(self, collection, data_list, is_show_log="true"):
        """
        插入多条记录，data_list样例格式[{"1":"2"},{"2":"3"},{"3": "4"}]
        这里只是抽象出，进行大规模的mongo插表,这里兼容了插入一条记录与多条记录
        :param collection:
        :param data_list:
        :param is_show_log: 用于设定是否要显示出插入成功的日志信息, 为false的时候就不显示成功的日志信息
        :return:
        """
        start = time.time()
        if not data_list:
            raise Exception("列表为空，不能进行后面的插表操作")
        table_id = None
        con = self.db[collection]
        record_num = len(data_list)
        try:
            if record_num > 5000000:
                for i in range(0, record_num, 4000000):
                    temp = data_list[i: i + 4000000]
                    con.insert_many(temp)
            else:
                if record_num >= 2:
                    con.insert_many(data_list)
                else:
                    table_id = con.insert_one(data_list[0]).inserted_id
        except Exception as e:
            if record_num >= 2:
                raise Exception("往{}插入多条记录失败{}".format(collection, e))
            else:
                raise Exception("往{}插入一条记录失败{}".format(collection, e))
        else:
            if is_show_log == 'true':
                if record_num >= 2:
                    print "往{}插入多条记录成功".format(collection)
                    # self.bind_object.logger.info("往{}插入多条记录成功".format(collection))
                else:
                    print "往{}插入一条记录成功".format(collection)
                    # self.bind_object.logger.info("往{}插入一条记录成功".format(collection))
                end = time.time()
                # self.bind_object.logger.info("文件导入mongo中花费的时间:{}".format(end - start))
                print ("文件导入mongo中花费的时间:{}".format(end - start))
            return table_id
    
    def check_objectid(self, id_):
        """
        用于检查并转成成ObjectID
        :param id_:
        :return:
        """
        if not isinstance(id_, ObjectId):
            if isinstance(id_, StringTypes):
                id_ = ObjectId(id_)
            else:
                raise Exception("id必须为ObjectId对象或其对应的字符串!")
        return id_