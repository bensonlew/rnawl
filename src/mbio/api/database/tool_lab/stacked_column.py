# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import datetime
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.ref_rna_v2.api_base import ApiBase
from bson.objectid import ObjectId
import types
import os
import re
from types import StringTypes
from collections import OrderedDict


class StackedColumn(ApiBase):
    def __init__(self, bind_object):
        super(StackedColumn, self).__init__(bind_object)
        self._project_type = 'tool_lab'


    def add_stacked_detail(self, file_path, new_otu_id=None, add_Algorithm=''):
        # if from_otu_id != 0 and not isinstance(from_otu_id, ObjectId):
        #     if isinstance(from_otu_id, StringTypes):
        #         from_otu_id = ObjectId(from_otu_id)
        #     else:
        #         self.bind_object.set_error("from_otu_table必须为ObjectId对象或其对应的字符串!", code="51000403")
        self.bind_object.logger.info("开始导入otu_detail表")
        if new_otu_id is None:
            name = "only_for_test"
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            # if type(params) == dict:
            #     params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_info = dict(
                # project_sn=project_sn,
                task_id="stacked",
                name=name,
                # exp_level=exp_level,
                # exp_type=exp_type.upper(),
                # method=quant_method,
                desc='stacked main table',
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                # params=params,
                version="v1",
                # sample_order=sample_order + group_order,
                status="start"
            )
            new_otu_id = self.create_db_table('stacked', [main_info])

        find_otu = self.db['stacked'].find_one({"_id": ObjectId(new_otu_id)})
        if find_otu:
            self.task_id = find_otu['task_id']
        else:
            self.bind_object.set_error("OTU_ID没有找到相关的主表信息", code="51000404")
        insert_data = list()
        cate=[]
        a=pd.read_table(file_path)
        b=a.T
        c = list(b.iloc[0])
        d = b.drop(b.index[0])
        e = d.to_dict("r")
        sample=list(d.index)
        # for p,i in enumerate(e):
        #     for n,t in enumerate(c):
        #         i[t]=i.pop(n)
        #     i["sample_id"]= sample[p]
        #     i["type"]="column"
        #     i["category"] = c,
        #     i['stacked_id'] = ObjectId(new_otu_id)
        insert_data=[]
        for p,i in enumerate(e):
            for n, t in enumerate(c):
                single_detail=dict()
                single_detail["name"] = sample[p]
                single_detail["type"] = "column"
                single_detail['stacked_id'] = ObjectId(new_otu_id)
                single_detail["category"] = t
                single_detail["value"] = i[n]
                insert_data.append(single_detail)
        # with open(file_path, 'rb') as r:
            # head = r.next().strip('\r\n')
            # head = re.split('\t', head)
            # new_head = head[1:]
            # for line in r:
            #     line = line.rstrip("\r\n")
            #     line = re.split('\t', line)
            #     sample_num = line[1:]
            #     classify_list = line[0]
            #     cate.append(classify_list)
            #     otu_detail = dict()
            #     otu_detail['stacked_id'] = ObjectId(new_otu_id)
            #     otu_detail["cate"] = classify_list
            #     for i in range(0, len(sample_num)):
            #         otu_detail[new_head[i]] = sample_num[i]
            #     otu_detail['task_id'] = self.task_id
            #     insert_data.append(otu_detail)
        try:
            collection = self.db['stacked_detail']
            collection.insert_many(insert_data)
            a=pd.read_table(file_path,index_col=0)
            lista=list(a.columns)
            dict_a={"name":"name","data":"value","category":"category","condition":{"type":"column"}}
            dict_b=json.dumps(dict_a, sort_keys=True, separators=(',', ':'))
            dict_c={"data":lista}
            dict_d = json.dumps(dict_c, sort_keys=True, separators=(',', ':'))
            # dict_b={"cate":cate}
            main_collection = self.db['stacked']  #guanqing.zou
            # main_collection.update({"_id":ObjectId(new_otu_id)}})  #guanqing 20180411
            main_collection.update({"_id": ObjectId(new_otu_id)}, {"$set": {"column_data":dict_b,"names":dict_d}})
        except Exception as e:
            self.bind_object.logger.error("导入stacked_detail表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入stacked_detail表格成功")

    def run1(self):
        self.add_stacked_detail("/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/test0409/taxa.percents.table.xls")


class TestFunction(unittest.TestCase):
    """
    测试导表函数
    """

    def test_mongo(test):
        from mbio.workflows.tool_lab.toollab_test_api import ToollabTestApiWorkflow
        from biocluster.wsheet import Sheet
        import random

        data = {
            # "id": "denovo_rna_v2" + str(random.randint(1,10000)),
            "id": "stacked",
            "project_sn": "stacked",
            "type": "workflow",
            "name": "tool_lab.toollab_test_api",
            "options": {
            },
        }
        wsheet = Sheet(data=data)
        wf = ToollabTestApiWorkflow(wsheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_END = False
        wf.test_api = wf.api.api("tool_lab.stacked_column")
        wf.test_api.run1()


if __name__ == '__main__':
    unittest.main()