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

project_type = "basic_molecular"
client = Config().get_mongo_client(mtype=project_type)
db = client[Config().get_mongo_dbname(project_type)]

class Lmplot(Base):
    def __init__(self, bind_object):
        super(Lmplot, self).__init__(bind_object)
        self._project_type = 'basic_molecular'
        
    def update_lmplot_record(self, main_id, target_path, s3_upload_path):
        txt_file = os.path.join(target_path)
        png_file = os.path.join(s3_upload_path)
        filelist = os.listdir(txt_file)
        lmplot_id = self.check_objectid(main_id)
        result = self.db["sg_lmplot"].find_one({"_id": lmplot_id})
        data_list = json.loads(result["format_data"])
        data_update = []
        for f in filelist:
            f = os.path.join(txt_file, f)
            dirname = os.path.dirname(f)
            baseName = os.path.basename(f)
            if dirname.endswith(os.sep):
                txt_name = dirname+baseName
            else:
                txt_name = dirname+os.sep+baseName
            with open(txt_name, 'r') as r:
                result_txt = r.readlines()
                slope = result_txt[0].strip().split('\t')[1]
                y_inter = result_txt[1].strip().split('\t')[1]
                efficiency = result_txt[2].strip().split('\t')[1]
                r2 = result_txt[3].strip().split('\t')[1]
            png = s3_upload_path + "{}_{}_line.png".format(baseName.split('_')[0],baseName.split('_')[1])
            input_txt_name = "{}_{}".format(baseName.split('_')[0],baseName.split('_')[1])
            input_png_name = "{}_{}".format(baseName.split('_')[0],baseName.split('_')[1])
            print input_png_name
            print input_txt_name
            input_slope = {input_txt_name :slope}
            input_y_inter = {input_txt_name :y_inter}
            input_efficiency = {input_txt_name :efficiency}
            input_r2 = {input_txt_name :r2}
            input_png = {input_txt_name :png}
            data_sort = list()
            data_dict = {}
            for data in data_list:
                sort_key = "{}_{}".format(data["format_data"][0]["gene"],data["id"])
                data_sort.append(sort_key)
                data_dict[sort_key] = data

            for in_key in input_slope.keys():
                data_dict[in_key]["slope"] = input_slope[in_key]

            for in_key in input_y_inter.keys():
                data_dict[in_key]["y_inter"] = input_y_inter[in_key]

            for in_key in input_r2.keys():
                data_dict[in_key]["r2"] = input_r2[in_key]

            for in_key in input_efficiency.keys():
                data_dict[in_key]["efficiency"] = input_efficiency[in_key]

            for in_key in input_png.keys():
                data_dict[in_key]["line_out_path"] = input_png[in_key]
            json_back = []
            for st_key in data_sort:
                json_back.append(data_dict[st_key])
            print(json.dumps(json_back))
            data_update.append(f)

            update_dict = {
                "lm_list" : json.dumps(json_back)
            }
            self.db['sg_lmplot'].update({"_id":self.check_objectid(main_id)},{"$set":update_dict}, upsert =True, multi= True)
        
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
