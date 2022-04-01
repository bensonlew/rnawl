# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from mbio.api.database.ref_rna_v2.api_base import ApiBase
from biocluster.api.database.base import report_check
import os
from bson.objectid import ObjectId
import unittest
import json
import pandas as pd

class RmatsModel(ApiBase):
    def __init__(self, bind_object):
        super(RmatsModel, self).__init__(bind_object)

    @report_check
    def add_rmats_model(self, output_dir, s3_output, main_id):
        insert_dict = {
            'pdf_files': [i for i in os.listdir(output_dir) if i.endswith('.pdf')],
            'png_files': [i for i in os.listdir(output_dir) if i.endswith('.png')],
            'graph_dir': s3_output,
            'main_id': ObjectId(main_id),
            'status': 'end'
        }
        self.update_db_record('sg_splicing_rmats_model', main_id, insert_dict=insert_dict)

    def add_rmats_model_sef(self, output_dir, main_id,group_dict):
        sample2group={}
        for group in group_dict:
            grouplist = group_dict[group]
            for sample in grouplist:
                sample2group[sample] = group
        insert_data=list()
        rmats_result_path = os.path.join(output_dir,"event_file.txt")
        plot_event_path = os.path.join(output_dir,"event_details.txt")
        rmat_df = pd.read_table(rmats_result_path,index_col=0)
        event_dict = rmat_df.to_dict("index")
        event_plot_dict=dict()
        event2path = {}
        with open(plot_event_path) as e:
            for line in e.readlines():
                event_id = line.strip().split(";")[0].split("\t")[1]
                event_title=line.strip().split(";")[1].split("\t")[1]
                event_file_path = line.strip().split(";")[2].split("\t")[1]
                event_plot_dict[event_id]=event_title
                event2path[event_id] = event_file_path
                data = {
                    'event_id': event_id,
                    'event_title' : event_title,
                    'Genomic': 'Genomic coordinate({})'.format(event_dict[event_id]["chr"]),
                    'strand': event_dict[event_id]["strand"],
                    'splicing_rmats_model_id': ObjectId(main_id),
                    'data_category': "title"
                }
                insert_data.append(data)
        for event_id in event_dict:
            detail_path = os.path.join(event2path[event_id],"detail")
            detail_files = os.listdir(detail_path)
            for detail_file in detail_files:
                sample_name = os.path.basename(detail_file).split("detail.txt")[0]
                sample_group = sample2group[sample_name]
                with open(os.path.join(detail_path,detail_file), 'r')as sd:
                    sample_detail = json.load(sd)

                connect_lines = []
                #这一段代码纯粹为了画图，标点，跨度更大的在下方展示
                connect_detail = sample_detail["line_connet"]
                if connect_detail:
                    jxn_range={}
                    jxn_up_down = {}
                    for jxn in sample_detail["line_connet"]:
                        left_point, right_point = jxn.split(":")
                        length = abs(int(right_point) - int(left_point))
                        jxn_range[jxn] = length
                    max_range = max([jxn_range[i] for i in jxn_range])
                    for jxn in sample_detail["line_connet"]:
                        left_point,right_point = jxn.split(":")
                        length = abs(int(right_point) - int(left_point))
                        left_hight= sample_detail["wiggle"][sample_detail["coord"].index(int(left_point))]
                        right_hight = sample_detail["wiggle"][sample_detail["coord"].index(int(right_point))]
                        reads_num =sample_detail["line_connet"][jxn]
                        if length < max_range:
                            connect_info=[[left_point,left_hight],[right_point,right_hight],reads_num]
                        else:
                            connect_info = [[left_point, 0], [right_point, 0], reads_num]
                        connect_lines.append(connect_info)
                data={
                    'event_id' : event_id,
                    'sample_name': sample_name,
                    'group' : sample2group[sample_name],
                    'sample_label' : sample_detail["sample_label"],
                    'wiggle': sample_detail["wiggle"],
                    'coord' : sample_detail["coord"],
                    'connect_line': connect_lines,
                    'data_category' : "sample",
                    'splicing_rmats_model_id' : ObjectId(main_id),
                }
                insert_data.append(data)
        try:
            self.create_db_table('sg_splicing_rmats_model_detail', insert_data)
            # collection = self.db['sg_splicing_rmats_model_detail']
            # collection.insert_many(insert_data)
        except Exception as e:
            self.bind_object.logger.info("导入splicing_rmats_module信息：%s出错:%s" % (plot_event_path, e))
        else:
            self.bind_object.logger.info("导入splicing_rmats_module信息：%s成功!" % (plot_event_path))

        # insert_dict = {
        #     'pdf_files': [i for i in os.listdir(output_dir) if i.endswith('.pdf')],
        #     'png_files': [i for i in os.listdir(output_dir) if i.endswith('.png')],
        #     'graph_dir': s3_output,
        #     'main_id': ObjectId(main_id),
        #     'status': 'end'
        # }
        event_list = [i for i in event_dict]
        update_dict = {
            # 'pdf_files': [i for i in os.listdir(output_dir) if i.endswith('.pdf')],
            # 'png_files': [i for i in os.listdir(output_dir) if i.endswith('.png')],
            # 'graph_dir': s3_output,
            'main_id': ObjectId(main_id),
            'status': 'end',
            "events":event_list,
            "dynamic" :"yes"
        }
        self.update_db_record('sg_splicing_rmats_model', main_id, insert_dict=update_dict)