# !/usr/bin/python
# -*- coding: utf-8 -*-
from bson.objectid import ObjectId
from api_base import ApiBase
import json
import pandas as pd
import datetime
import glob
import os


class TfbsPredict(ApiBase):
    def __init__(self, bind_object):
        super(TfbsPredict, self).__init__(bind_object)

    def add_tfbs_predict_detail(self, output_dir, main_id):
        result_file = glob.glob(output_dir + '/tf*.xls')
        # get gene_detail
        gene_info = os.path.join(os.path.dirname(os.path.dirname(output_dir)), 'seq_annot.xls')
        gene_info_pd = pd.read_table(gene_info, header=0,)
        #gene_info_pd = gene_info_pd.loc[:, ["gene_id", "gene_name", "gene_desc"]].drop_duplicates()
        gene_info_pd = gene_info_pd.loc[:, ["gene_id", "gene_name"]].drop_duplicates()
        gene_info_pd.set_index("gene_id", inplace=True)
        for each in result_file:
            result_pd = pd.read_table(each, header=0)
            if "tf_geneid" in result_pd.columns:
                result_pd = result_pd.join(gene_info_pd, on="tf_geneid")
            if "target_id" in result_pd.columns:
                result_pd = result_pd.join(gene_info_pd, on="target_id", lsuffix="_tf", rsuffix="_target")
            result_pd = result_pd.fillna('-')
            result_pd['tfbs_predict_id'] = ObjectId(main_id)
            result_dict_list = result_pd.to_dict("records")
            if each.endswith("_predict.xls"):
                self.create_db_table('tfbs_predict_detail', result_dict_list)
            elif each.endswith("tf_target_number.xls"):
                self.create_db_table('tfbs_predict_stat_detail', result_dict_list)
        self.update_db_record('tfbs_predict', main_id, status="end", )


if __name__ == '__main__':
    TfbsPredict = TfbsPredict(None)
    TfbsPredict.add_tfbs_predict_detail("/mnt/ilustre..../WgcnaModule", "5aaf7efca4e1af606ad5e527")
