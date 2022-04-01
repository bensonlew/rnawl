# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'xueqinwen'

from collections import OrderedDict
from bson import SON
from biocluster.config import Config
from io import StringIO
import types
import re
import datetime
import os
import json
from api_base import ApiBase

class BacIdentity(ApiBase):
    """
    这是菌鉴流程在测试机上测试的导表
    对应的是测试的sanger_tool_lab的mongo库
    """
    def __init__(self, bind_object):
        super(BacIdentity, self).__init__(bind_object)
        self._db_name = 'sanger_tool_lab'

    def update_bac_db_record(self, main_id, s3_upload_dir,final_dir):
        """
        更新主表字段
        """
        update_dict = {
            "s3": s3_upload_dir,
            "final_name":final_dir
        }
        self.update_db_record(
            "bac_identity",{"_id":self.check_objectid(main_id)}, update_dict
        )

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


    