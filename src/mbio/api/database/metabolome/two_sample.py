# -*- coding: utf-8 -*-


from biocluster.api.database.base import Base, report_check
import os
import json
from collections import defaultdict
import datetime
from bson.objectid import ObjectId
from types import StringTypes
import pandas as pd
import glob

class TwoSample(Base):
    def __init__(self, bind_object):
        super(TwoSample, self).__init__(bind_object)
        self._project_type = "metabolome"

    @report_check
    def add_two_sample(self, name=None, params=None, metab_table_id=None):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "desc": "运行结束",
            "name": name,
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "metab_table_id" : metab_table_id,
            "version" :'3.0'

        }
        collection = self.db["diff_sample"]
        inserted_id = collection.insert_one(insert_data).inserted_id
        collection.update_one({'_id': inserted_id}, {'$set': {'main_id': inserted_id}})
        return inserted_id

    @report_check
    def add_two_sample_detail(self, result_path, main_id, table_type=None):
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)
        data_list = []
        diff_groups = os.listdir(result_path)
        for diff_group in diff_groups:
            result_files = glob.glob(result_path+'/'+diff_group+'/*_result.xls')
            for file in result_files:
                data = pd.read_table(file,sep='\t',index_col=0)
                for index in data.index:
                    insert_data = {
                        "diff_id" : main_id,
                        "metab_id" : index,
                        "fdr" : data.loc[index]['corrected_pvalue'],
                        "fc" : data.loc[index]['odds_ratio'],
                        "p_value" : data.loc[index]['pvalue'],
                        "diff_group" : diff_group,
                        "metab" :data.loc[index]['Metabolite']
                    }
                    data_list.append(insert_data)
        if table_type:
            for data in data_list:
                data['table_type'] = table_type

        try:
            collection = self.db["diff_sample_detail"]
            collection.insert_many(data_list)
            self.db['diff_sample'].update({"_id":main_id},{"$set": {"main_id": main_id}})
        except Exception as e:
            self.bind_object.logger.info("导入diff_sample_detail数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入diff_sample_detail数据成功")

