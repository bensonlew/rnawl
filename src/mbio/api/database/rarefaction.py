# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
# last modified guhaidong 20171116
from biocluster.api.database.base import Base, report_check
import os
import datetime
from bson.objectid import ObjectId
from types import StringTypes
import json
# from biocluster.config import Config
from mainapp.libs.param_pack import group_detail_sort
import copy


class Rarefaction(Base):
    def __init__(self, bind_object):
        super(Rarefaction, self).__init__(bind_object)
        self._project_type = 'meta'
        # self._db_name = Config().MONGODB
        self.category_x = []
        self.category_y = []

    @report_check
    def add_rarefaction_detail(self, rare_id, file_path):
        if not isinstance(rare_id, ObjectId):
            if isinstance(rare_id, StringTypes):
                rare_id = ObjectId(rare_id)
            else:
                self.bind_object.set_error("rarefaction_id必须为ObjectId对象或其对应的字符串!", code="51005501")
        collection_first = self.db['sg_alpha_rarefaction_curve']
        result = collection_first.find_one({"_id": rare_id})  # readd by guhaidong 20171116
        task_id = result['task_id']  # readd by guhaidong 20171116
        rare_paths = os.listdir(file_path)
        rare_detail = []
        category_x = []
        max_counts = []
        for rare_path in rare_paths:
            # print os.path.join(file_path,rare_path)
            files = os.listdir(os.path.join(file_path, rare_path))
            # print files
            fs_path = []
            for f in files:
                fs_path.append(os.path.join(file_path, rare_path, f))
                # self.bind_object.logger.error(fs_path)
            for fs in fs_path:
                rarefaction = []
                sample_name = fs.split('.')[1]
                with open(fs) as f:
                    while True:
                        line = f.readline().strip('\n')
                        if not line:
                            break
                        line_data = line.split("\t")
                        # print line_data
                        category_x.append(line_data[0])
                        self.category_y.append(float(line_data[1]))
                        my_dic = dict()
                        my_dic["column"] = line_data[0]
                        my_dic["value"] = line_data[1]
                        rarefaction.append(my_dic)
                    rarefaction.pop(0)
                    max_counts.append(len(rarefaction))
                    # print rarefaction
                    insert_data = {
                        "rarefaction_curve_id": rare_id,
                        # "task_id": task_id,
                        "index_type": rare_path,
                        "specimen_name": sample_name,
                        "json_value": rarefaction
                    }
                    # print insert_data
                    rare_detail.append(insert_data)
        x = copy.copy(category_x)
        for i in x:
            if i == 'numsampled':
                category_x.remove(i)
            else:
                self.category_x.append(int(i))
        try:
            collection = self.db['sg_alpha_rarefaction_curve_detail']
            collection.insert_many(rare_detail)
            # collection_first = self.db['sg_alpha_rarefaction_curve']
            collection_first.update({"_id": ObjectId(rare_id), "task_id": task_id},
                                    {"$set": {"max_count": max(max_counts), "category_x": max(self.category_x), "status": "end"}})
            #collection_first.update({"_id": ObjectId(rare_id)},{"$set": {"main_id":  ObjectId(rare_id)}})
            # add task_id by guhaidong 20171116
        except Exception as e:
            self.bind_object.logger.error("导入rare_detail表格{}信息出错:{}".format(file_path, e))
            self.bind_object.set_error("导入rare_detail表格出错", code="51005502")
        else:
            self.bind_object.logger.info("导入rare_detail表格{}成功".format(file_path))
        # return max(self.category_x)

    @report_check
    def add_rare_table(self, file_path, level, otu_id=None, task_id=None, name=None, params=None, spname_spid=None):
        if level not in range(1, 10):
            self.bind_object.logger.error("level参数%s为不在允许范围内!" % level)
            self.bind_object.set_error("level参数不在允许的范围内", code="51005503")
        if task_id is None:
            task_id = self.bind_object.sheet.id
        if otu_id:
            if not isinstance(otu_id, ObjectId):
                otu_id = ObjectId(otu_id)
            params['otu_id'] = str(otu_id)  # otu_id在再metabase中不可用
        if spname_spid:
            group_detail = {'All': [str(i) for i in spname_spid.values()]}
            params['group_detail'] = group_detail_sort(group_detail)
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": task_id,
            "otu_id": otu_id,
            "name": name if name else "Rarefaction_Origin",
            "level_id": level,
            "status": "start",
            "desc": "",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':')),
            # "group_id": group_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        # if params is not None:
            # insert_data['params'] = json.dumps(params, sort_keys=True, separators=(',', ':'))
        collection = self.db["sg_alpha_rarefaction_curve"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        self.add_rarefaction_detail(inserted_id, file_path)
        return inserted_id
