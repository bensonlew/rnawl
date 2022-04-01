# -*- coding: utf-8 -*-



import os
from bson.objectid import ObjectId
import pandas as pd
import glob
from api_base import ApiBase
import json

class Fisher(ApiBase):
    def __init__(self, bind_object):
        super(Fisher, self).__init__(bind_object)

    def add_two_sample_detail(self, result_path, main_id,sample1, sample2):
        if not isinstance(main_id,ObjectId):
            main_id = ObjectId(main_id)
        data_list = []
        sig_list =[]
        pvalue_text_list = []
        result_files = glob.glob(result_path+'/*_result.xls')
        file = result_files[0]
        file_ci = glob.glob(result_path+'/*_CI.xls')[0]
        data = pd.read_table(file,sep='\t',index_col=0)

        data_ci = pd.read_table(file_ci, sep='\t',index_col=0)
        data = data.join(data_ci)
        data.fillna(0,inplace=True)
        data['p_rank'] = data['pvalue'].rank(ascending=True,method='first')

        id = 0
        for index in data.index:
            id +=1
            tmp = data.loc[index]
            fc = round(tmp['odds_ratio'],5)
            if  tmp['odds_ratio']==0:
                d_fc = 0
            else:
                d_fc = round(float(1)/tmp['odds_ratio'],5)
            mix_fc = max([fc,d_fc])

            insert_data = {
                "diff_id" : main_id,
                "name" : index,
                "fc" : fc,
                "d_fc" : d_fc,
                "mix_fc": mix_fc,
                "rank" : id,
                "p_rank" : tmp['p_rank'],
                sample1 : tmp[sample1],
                sample2 : tmp[sample2],
                "p_value" : tmp['pvalue'],
                "q_value" : tmp['corrected_pvalue'],
                "effect" : tmp['effectsize'],
                "low_ci"  : tmp["lowerCI"],
                "up_ci"  : tmp["upperCI"]
            }
            data_list.append(insert_data)

            sig_data = {
                "diff_id" : main_id,
                "name": index,
                "value" :  tmp['pvalue'],
                "x" : index, #"{0}|{1}".format(sample1, sample2),
                "type": "significant",
                "factor": "0.05|0.01|0.001"
            }
            sig_list.append(sig_data)

            text_data = {
                "diff_id" : main_id,
                "name": index,
                "x" : index,
                "text" : tmp['pvalue'],
                "type" : "text"
            }
            pvalue_text_list.append(text_data)

        try:
            collection = self.db["fisher_detail"]
            collection.insert_many(data_list)
            self.db["fisher_sig"].insert_many(sig_list)
            self.db["fisher_sig"].insert_many(pvalue_text_list)

            table_data = {"column":[{"filter": "false", "field": "name", "type": "string", "sort": "false", "title": "Name"},
                                    {"filter": "false", "field": "p_value", "type": "string", "sort": "false", "title": "Pvalue"},
                                    {"filter": "false", "field": "q_value", "type": "string", "sort": "false", "title": "Corrected Pvalue"},
                                    {"filter": "false", "field": "fc", "type": "string", "sort": "false", "title": "Fold Change(%s/%s)"%(sample1,sample2)},
                                    {"filter": "false", "field": sample1, "type": "string", "sort": "false", "title": sample1},
                                    {"filter": "false", "field": sample2, "type": "string", "sort": "false", "title":sample2},
                                    {"filter": "false", "field": "low_ci", "type": "string", "sort": "false", "title": "Lower ci"},
                                    {"filter": "false", "field": "up_ci", "type": "string", "sort": "false", "title":"Upper ci"},
                                    {"filter": "false", "field": "effect", "type": "string", "sort": "false", "title":"Effect size"}],
                          "condition": {}}
            table_data = json.dumps(table_data)
            column_data = json.dumps({"name":"name","data":[sample1,sample2]})
            significant_data = json.dumps({"name": "name","value":"value", "condition":{"type": "significant"}})
            pvalue_text_data = json.dumps({"name": "name","condition":{"type":"text"}})
            self.db["fisher"].update({"_id": main_id}, {"$set": {
                "main_id":main_id,
                "table_data":table_data,
                "column_data":column_data,
                "significant_data":significant_data,
                "pvalue_text_data" : pvalue_text_data
            }})
        except Exception as e:
            self.bind_object.logger.info("导入明细数据出错:%s" % e)
        else:
            self.bind_object.logger.info("导入明细数据成功")

