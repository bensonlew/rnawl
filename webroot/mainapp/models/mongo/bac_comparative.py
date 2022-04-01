# -*- coding: utf-8 -*-
#__author__ = 'hao.gao'

import types
import json
import random, os, re
from bson import SON
from bson.objectid import ObjectId
from biocluster.config import Config
from mainapp.models.workflow import Workflow
from mainapp.models.mongo.meta import Meta
from biocluster.file import exists,list_dir


class BacComparative(Meta):
    """
    获取mongo字段
    """
    def __init__(self, bind_object=None):
        self._bind_object = bind_object
        super(BacComparative, self).__init__(self._bind_object)
        self._project_type = 'bac_comparative'

    def get_projectsn(self, task_id):
        """
        获取project_sn
        :param task_id:
        :return:
        """
        collection = self.db['sg_task']
        result = collection.find_one({"task_id": task_id})
        if not result:
            raise Exception('{}没有响应的id:{}!'.format(collection, task_id))
        project_sn = result["project_sn"]
        return project_sn

    def get_panfile(self, task_id, pan_id):
        """
        获取cluster_path
        :param task_id,pan_id
        :return:
        """
        collection = self.db['pan']
        result = collection.find_one({"task_id": task_id, "main_id": pan_id})
        if not result:
            raise Exception('{}没有响应的id:{},pan_id:{}!'.format(collection, task_id, pan_id))
        cluster_path = result["cluster_path"]
        return cluster_path

    def new_info_to_status(self, table_id, add_dict):
        '''
        sg_status 中增加进行字段 by xieshichang
        '''
        collection = self.db['sg_status']
        sg_status = collection.find_one({'table_id': table_id})
        if not sg_status:
            raise Exception('{}没有响应的id:{}!'.format(collection, 'sg_status'))
        collection.update({'table_id': table_id}, {'$set': add_dict})

    def get_value_of_key(self, table, table_id, mongo_key):
        '''
        根据表中字段的值
        '''
        info = self.get_main_info(table_id, table)
        return info[mongo_key]

    def get_pan_category(self, table_id):
        pan_group = self.get_main_info(table_id, 'pan_group')
        category = pan_group['category']
        value_type = pan_group['value_type']
        category_name = pan_group['category_names']
        # category_scheme = pan_group['category_scheme']
        if not value_type or (value_type in ["true", "True"]):
            # all_value = [v2 for k1, v1 in category.items() for k2, v2 in v1.items()]
            max_v = 100 ## 转成百分比内容
            category = {k: {"min": float(v['min'])/max_v, "max": float(v['max'])/max_v} for k, v in category.items()}
        else:
            category = {k: {"min": v['min'], "max":v['max']} for k, v in category.items()}

        return (json.dumps(category), value_type, category_name)

    def get_pan_category2(self, table_id, sample_num=None):
        pan_group = self.get_main_info(table_id, 'pan_group')
        category = pan_group['category']
        value_type = pan_group['value_type']
        all_value = []
        if not value_type or value_type in ['false', 'False']: # 20200508
            for k1, v1 in category.items():
                for k2, v2 in v1.items():
                    category[k1][k2] = eval(str(v2).replace('N', str(sample_num)))
                    all_value.append(category[k1][k2])
            #all_value = [float(v2) for k1, v1 in category.items() for k2, v2 in v1.items()]
            max_v = max(all_value) / 100.0
            c = [[k, float(v['min'])/max_v, float(v['max'])/max_v] for k, v in category.items()]
            c.sort(key=lambda x: x[1])
            for i in range(len(c) -1):
                c[i][2] = c[i+1][1]
            category = {a[0]: [a[1], a[2]] for a in c}
        else:
            category = {k: [float(v['min']), float(v['max'])] for k, v in category.items()}

        return json.dumps(category)

    def get_annofile(self, task_id):
        """
        获取anno_file总览表
        :param task_id
        :return:
        """
        collection = self.db['anno_summary']
        result = collection.find_one({"task_id": task_id, "name": "summary_origin"})
        if not result:
            raise Exception('{}没有响应的id:{}!'.format(collection, task_id))
        summary_path = result["summary_path"]
        return summary_path
