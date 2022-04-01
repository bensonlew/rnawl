# !/usr/bin/python
# -*- coding: utf-8 -*-
from bson.objectid import ObjectId
from api_base import ApiBase
import json
import pandas as pd
import datetime
import glob


class TfPredict(ApiBase):
    def __init__(self, bind_object):
        super(TfPredict, self).__init__(bind_object)

    def add_tf_predict_detail(self, output_dir, main_id):
        find_file = glob.glob(output_dir + '/final_tf_predict.xls')
        if find_file:
            result_file = find_file[0]
            result_pd = pd.read_table(result_file, header=0).drop_duplicates()
            result_pd = result_pd.fillna("")
            families = list(set(result_pd['family']))
            result_pd['tf_predict_id'] = ObjectId(main_id)
            result_dict_list = result_pd.to_dict("records")
            self.create_db_table('sg_tf_predict_detail', result_dict_list)
            self.update_db_record('sg_tf_predict', main_id, status="end", families=families)


if __name__ == '__main__':
    TfPredict = TfPredict(None)
    TfPredict.add_tf_predict_detail("/mnt/ilustre..../WgcnaModule", "5aaf7efca4e1af606ad5e527")
