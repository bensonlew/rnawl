# !/usr/bin/python
# -*- coding: utf-8 -*-
from bson.objectid import ObjectId
from api_base import ApiBase
import json
import pandas as pd
import datetime
import glob


class TfStat(ApiBase):
    def __init__(self, bind_object):
        super(TfStat, self).__init__(bind_object)

    def add_tf_stat_detail(self, output_dir, main_id):
        result_file = glob.glob(output_dir + '/*.txt')
        for each_file in result_file:
            result_pd = pd.read_table(each_file, header=0)
            if '_bar_' in each_file:
                # tmp_dict = dict()
                # tmp_dict['family'] = list(result_pd['family'])
                # tmp_dict['gene_num'] = list(result_pd['gene_num'])
                # tmp_dict['transcript_num'] = list(result_pd['transcript_num'])
                # tmp_dict['tf_stat_id'] = ObjectId(main_id)
                # self.create_db_table('tf_stat_bar_detail', [tmp_dict])
                result_pd['tf_stat_id'] = ObjectId(main_id)
                result_dict_list = result_pd.to_dict("records")
                self.create_db_table('tf_stat_bar_detail', result_dict_list)
            elif '_circos_' in each_file:
                result_pd['tf_stat_id'] = ObjectId(main_id)
                result_dict_list = result_pd.to_dict("records")
                self.create_db_table('tf_stat_circos_detail', result_dict_list)
        self.update_db_record('tf_stat', main_id, status="end", )


if __name__ == '__main__':
    TfStat = TfStat(None)
    TfStat.add_tf_stat_detail("/mnt/ilustre..../WgcnaModule", "5aaf7efca4e1af606ad5e527")
