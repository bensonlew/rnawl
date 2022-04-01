# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modify:20210924
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
from biocluster.config import Config
import datetime
import json
import os
import sys
from mbio.api.database.ref_rna_v2.api_base import ApiBase
from pymongo import ASCENDING, DESCENDING, HASHED

class EnsureIndex(ApiBase):
    def __init__(self, bind_object):
        super(EnsureIndex, self).__init__(bind_object)


    def get_index(self):
        '''
        获取所有collection 对应的index， 最好在线上运行
        '''
        self.index_dict = dict()
        cols = self.db.list_collection_names()
        for col in cols:
            index = self.db[col].index_information()
            index.pop("_id_")
            self.index_dict[col] = index
        with open("all_index.json", "w") as f:
            f.write(json.dumps(self.index_dict, indent=4, ensure_ascii=False).encode('utf8').decode())

    def ensure_all_indexes(self, json_file=None):
        script_path = os.path.split(os.path.realpath(__file__))[0]
        if not json_file:
            index_json_file = os.path.join(script_path, "all_index.json")

        with open(index_json_file, "r") as f:
            index_dict = json.load(f)

        for col, indexes in index_dict.items():
            for name, index in indexes.items():
                self.create_index(col, name, index)

    def create_index(self, collection, name, index):
        if index:
            col = self.db[collection]
            index_create = []
            for k in index["key"]:
                if k[1] == 1:
                    index_create.append((k[0], ASCENDING))
                elif k[1] == -1:
                    index_create.append((k[0], DESCENDING))
                elif k[1] == "hashed":
                    index_create.append((k[0], HASHED))
                else:
                    print("unknown index {} {}".format(collection, name))
            col.create_index(index_create, name=name)
            print("{}的{}索引{}构建完成！".format(collection, name, index))
            #self.bind_object.logger.info("{}的{}索引{}构建完成！".format(collection, name, index))
        else:
            print("传入的index值为空，不进行建索引！")
            # self.bind_object.logger.info("传入的index值为空，不进行建索引！")


if __name__ == '__main__':
    os.environ["current_mode"]="workflow"
    os.environ["NTM_PORT"]="7322"
    os.environ["WFM_PORT"]="7321"

    report = EnsureIndex(None)
    # report._config.DBVersion = 0
    # # report._config.ProjectID = "v0afleo6p1211n885l0oo0av7v"
    # report.get_report_model(task_id=sys.argv[1])

    if sys.argv[1] == 'get_index':
        # report._config.ProjectID = "v0afleo6p1211n885l0oo0av7v"
        report._config.DBVersion = 1
        report.get_index()
    elif sys.argv[1] == 'create_index':
        report._config.DBVersion = 1
        report.ensure_all_indexes()
        pass
