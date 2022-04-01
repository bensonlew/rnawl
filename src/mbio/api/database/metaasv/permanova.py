# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:2018.04.16
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
# from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId


class Permanova(Base):
    def __init__(self, bind_object):
        super(Permanova, self).__init__(bind_object)
        self._project_type = "metaasv"

    @report_check
    def add_permanova(self, dir_path, main=False, task_id=None, main_id=None,
                      anno_type=None, params=None, name=None):
        if main:
            if task_id is None:
                task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else 'PERMANOVA_Origin',
                'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
                'status': 'end',
                'anno_type': anno_type
            }
            collection = self.db['permanova']
            # 将主表名称写在这里
            main_id = collection.insert_one(insert_data).inserted_id
        else:
            if main_id is None:
                self.bind_object.set_error("main为False时需提供main_id!")
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
            _main_collection = self.db['permanova']
            result = _main_collection.find_one({'_id': main_id})
            if result:
                if os.path.exists(dir_path):
                    one_factor = os.path.join(dir_path, "One-factor_PERMANOVA.xls")
                    self.add_permanova_table(one_factor, main_id, "single")
                    multi_factor = os.path.join(dir_path, "Multi-factor_PERMANOVA.xls")
                    self.add_permanova_table(multi_factor, main_id, "multiple")
            else:
                self.bind_object.logger.error('提供的_id：%s在分析记录中无法找到表, taskid: %s' % (str(main_id), self.task_id))
                self.bind_object.set_error("找不到表")
            return main_id

    def add_permanova_table(self,file_path, main_id, type):
        """
        导入详情表permanova_table
        :param file_path:
        :param main_id:
        :param type:
        :return:
        """

        data_list = []
        with open(file_path, 'r') as f:
            lines = f.readlines()
            collection_detail = self.db["permanova_table"]
            for line in lines[1:]:
                if type in ["multiple"]:
                    data = [("permanova_id", main_id)]
                    line = line.strip().split('\t')
                    data.extend([("name", line[0]),
                                 ("type", "multiple"),
                                 ("df", line[1]),
                                 ("sum_of_sqs", line[2]),
                                 ("meanseq", line[3]),
                                 ("f_model", line[4]),
                                 ("r2", line[5]),
                                 ("p_value", line[6]),
                                 ("p_adjust", line[7])])
                    data_son = SON(data)
                    data_list.append(data_son)
                elif type in ["single"]:
                    data = [("permanova_id", main_id)]
                    line = line.strip().split('\t')
                    data.extend([("name", line[0]),
                                 ("type", "single"),
                                 ("sum_of_sqs", line[1]),
                                 ("meanseq", line[2]),
                                 ("f_model", line[3]),
                                 ("r2", line[4]),
                                 ("p_value", line[5]),
                                 ("p_adjust", line[6])])
                    data_son = SON(data)
                    data_list.append(data_son)
        try:
            collection_detail.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入permanova%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入permanova信息出错")
        else:
            self.bind_object.logger.info("导入permanova%s信息成功!" % file_path)
        try:
            main_collection = self.db["permanova"]
            settled_params = {"software" : "R-3.3.1 (vegan)"}
            settled_params_json = json.dumps(settled_params, sort_keys=True, separators=(',', ':'))
            if type in ["multiple"]:
                multiple_table_data = {
                        "table_data": ["name","df", "sum_of_sqs", "meanseq", "f_model", "r2","p_value"],
                        "condition": {"type": "multiple"}
                    }
                table_data_json = json.dumps(multiple_table_data, sort_keys=True, separators=(',', ':'))
                main_collection.update({"_id":main_id},{"$set":{"main_id": main_id,
                                                                "settled_params":settled_params_json,
                                                                "multiple_table_data": table_data_json}})
            elif type in ["single"]:
                single_table_data = {
                        "table_data": ["name", "sum_of_sqs", "meanseq", "f_model", "r2","p_value"],
                        "condition": {"type": "single"}
                    }
                table_data_json = json.dumps(single_table_data, sort_keys=True, separators=(',', ':'))
                main_collection.update({"_id":main_id},{"$set":{"main_id": main_id,
                                                                "settled_params":settled_params_json,
                                                                "single_table_data": table_data_json}})
        except Exception as e:
            self.bind_object.logger.error("导入permanova%s信息出错:%s" % (file_path, e))
            self.bind_object.set_error("导入permanova信息出错")
        else:
            self.bind_object.logger.info("导入permanova%s信息成功!" % file_path)
        return main_id