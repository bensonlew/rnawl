# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'
# last_modify:20210202
from biocluster.api.database.base import Base, report_check
from bson.objectid import ObjectId
from types import StringTypes
from bson.son import SON
from biocluster.config import Config
import datetime
import json
import os
import pandas as pd
import glob
import re
from mbio.api.database.ref_rna_v2.api_base import ApiBase

class DiffRadar(ApiBase):
    def __init__(self, bind_object):
        super(DiffRadar, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_radar_result(self, main_id, radar_file):
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        text_dict_raw = {"name": "name", "condition": {"type": "text"},"radar_type":"radar_type","group":"category"}
        scatter_data_raw = {"name":"name","data":["seq_id","log2fc"],"category":"regulate","condition":{"type":"scatter"},"size":"size"}
        # line_data_raw = {"name":"name", "category":"category","condition":{"type":"c_peak"}}
        line_data_raw = {"name": "name", "category": "category", "condition": {"type": "line"}}

        text_dict = json.dumps(text_dict_raw, sort_keys=True, separators=(',', ':'))
        scatter_data = json.dumps(scatter_data_raw, sort_keys=True, separators=(',', ':'))
        line_data = json.dumps(line_data_raw, sort_keys=True, separators=(',', ':'))




        #接下来开始准备数据
        df = pd.read_table(radar_file, header=0, sep='\t')
        # df = df.round(2)
        for i in [1,2,3]:
            df.iloc[:, i] = df.iloc[:, i].apply(lambda x: round(float(x + 0.0001), 2))

        #准备真正的外圈信息
        text0_df = df.iloc[:, [0]]
        text0_df.columns = ["name"]
        text0_df["text"] =  text0_df["name"]
        text0_df["radar_type"] = "outset"
        text0_df["type"] = "text"
        text0_df["diff_radar_id"] = main_id
        text0_detail = text0_df.to_dict("r")

        #首先准备外圈文字数据信息
        text1_df = df.iloc[:,[0,3]]
        text1_df.columns = ["name","text"]
        text1_df["radar_type"] = "outer"
        text1_df["type"] = "text"
        text1_df["diff_radar_id"] = main_id
        text1_detail = text1_df.to_dict("r")

        # 再准备散点信息
        scatter_df = df.iloc[:,[0,3]]
        scatter_df.columns = ["seq_id", "log2fc"]
        scatter_df["regulate"] = scatter_df["log2fc"].apply(lambda x :"up" if x >0 else "down")
        scatter_df["size"] = abs(scatter_df["log2fc"])
        # scatter_df["name"] = scatter_df.index
        scatter_df["name"] = scatter_df["seq_id"]
        scatter_df["type"] = "scatter"
        scatter_df["diff_radar_id"] = main_id
        scatter_detail = scatter_df.to_dict("r")

        # 再准备内圈文字信息
        control_group_name =  df.columns[1]
        compare_group_name = df.columns[2]
        text2_df = df.iloc[:,]
        text2_df.columns = ["x",control_group_name,compare_group_name,"y"]
        text2_df_split_1 = text2_df[["x",control_group_name]]
        text2_df_split_1.columns = ["x", "y"]
        text2_df_split_1["category"] = "abundance of {}".format(control_group_name)
        text2_df_split_2 = text2_df[["x", compare_group_name]]
        text2_df_split_2.columns = ["x", "y"]
        text2_df_split_2["category"] = "abundance of {}".format(compare_group_name)
        text2_final_df = pd.concat([text2_df_split_1,text2_df_split_2])
        text2_final_df["name"] = text2_final_df.index
        text2_final_df["type"] = "text"
        text2_final_df["diff_radar_id"] = main_id
        text2_final_df["radar_type"] = "inner"
        text2_detail = text2_final_df.to_dict("r")

        #最后准备内部线型图
        control_group_name = df.columns[1]
        compare_group_name = df.columns[2]
        peak_df = df.iloc[:, ]
        peak_df.columns = ["x", control_group_name, compare_group_name, "y"]
        peak_df_split_1 = peak_df[["x", control_group_name]]
        peak_df_split_1.columns = ["x", "y"]
        peak_df_split_1["category"] = "abundance of {}".format(control_group_name)
        peak_df_split_2 = peak_df[["x", compare_group_name]]
        peak_df_split_2.columns = ["x", "y"]
        peak_df_split_2["category"] = "abundance of {}".format(compare_group_name)
        peak_final_df = pd.concat([peak_df_split_1, peak_df_split_2])
        peak_final_df["name"] = text2_final_df.index
        # peak_final_df["type"] = "c_peak"
        peak_final_df["type"] = "line"
        peak_final_df["diff_radar_id"] = main_id
        peak_detail = peak_final_df.to_dict("r")
        self.create_db_table('diff_radar_detail', text0_detail)
        self.create_db_table('diff_radar_detail', text1_detail)
        self.create_db_table('diff_radar_detail', scatter_detail)
        self.create_db_table('diff_radar_detail', text2_detail)
        self.create_db_table('diff_radar_detail', peak_detail)

        scatter_group_dict = {"data":list(scatter_df["regulate"].drop_duplicates())}
        scatter_group_data = json.dumps(scatter_group_dict, sort_keys=True, separators=(',', ':'))
        self.update_db_record('diff_radar', main_id, main_id=main_id, text_data=text_dict,line_data = line_data,
                              scatter_data=scatter_data,scatter_group=scatter_group_data ,status = "end")




