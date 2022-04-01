# -*- coding: utf-8 -*-
# __author__ = 'xuxi'
# last_modify:20210729
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
from api_base import ApiBase
from mainapp.models.mongo.tool_lab import ToolLab

class Bicluster(ApiBase):
	def __init__(self, bind_object):
		super(Bicluster, self).__init__(bind_object)
		self._project_type = 'tool_lab'

	def add_bicluster(self, pdf_num, tool_output_dir, s3_dir, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
		if not main_id:
			name = "bicluster"+'_'
			time_now = datetime.datetime.now()
			name += time_now.strftime("%Y%m%d_%H%M%S")
			main_info = dict(
				project_sn=project_sn,
				task_id=task_id,
				version="v1",
				name=name,
				created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
				desc='bicluster',
				params= params if params else "null",
				status="start",
			)
			main_id = ToolLab().create_db_table('sg_bicluster', [main_info])
		else:
			main_id = ObjectId(main_id)

		pdf_name_list = ["sub_cluster_{}.pdf".format(str(i+1)) for i in range(pdf_num)]
		tsv_name_list = ["sub_cluster_{}.tsv".format(str(i+1)) for i in range(pdf_num)]
		pdf_name_list_existed = []
		pdf_dict = {}


		tsv_info = {'column': {}, 'condition': {}}
		for tsv_name in tsv_name_list:
			if os.path.exists(os.path.join(tool_output_dir ,tsv_name)):
				data_df = pd.read_table(os.path.join(tool_output_dir ,tsv_name), header=0)
				this_tsv_columns_ = list(data_df.columns)[1:]
				this_tsv_columns = [{"filter": "false","field": i,"title": i,"type": 'float',"sort": "false"} for i in this_tsv_columns_]
				this_tsv_columns.insert(0,{"filter": "false","field": 'gene',"title": 'gene',"type": 'string',"sort": "false"})
				tsv_info['column'][tsv_name.replace('.tsv','')] = this_tsv_columns

		for pdf_name in pdf_name_list:
			if os.path.exists(os.path.join(tool_output_dir ,pdf_name)):
				pdf_name_list_existed.append(pdf_name.replace('.pdf',''))
				pdf_dict[pdf_name.replace('.pdf','')] = os.path.join(s3_dir, pdf_name)
		try:
			update_dict = {
				"status":"end",
			 	"pdf_num":pdf_num,
			 	"pdf_list":json.dumps({"data":pdf_name_list_existed}, sort_keys=True, separators=(',', ':')),
			 	"pdf_dict":pdf_dict,
			 	"detail_data":json.dumps(tsv_info, sort_keys=True, separators=(',', ':'))
			 }
			self.update_db_record('sg_bicluster', {"_id": ObjectId(main_id)}, update_dict)
		except Exception as e:
			self.bind_object.logger.info("导入双聚类数据出错:%s" % e)
		else:
			self.bind_object.logger.info("导入双聚类数据成功")

		for tsv_name in tsv_name_list:
			if os.path.exists(os.path.join(tool_output_dir ,tsv_name)):
				self.add_tsv_detail(os.path.join(tool_output_dir ,tsv_name),bicluster_id=main_id, sub_cluster=tsv_name.replace('.tsv',''))

	def add_tsv_detail(self, tsv, bicluster_id, sub_cluster):
		data_df = pd.read_table(tsv, header=0)
		data_df["bicluster_id"] = bicluster_id
		data_df['sub_cluster'] = sub_cluster
		data_df.rename( columns={'Unnamed: 0':'gene'}, inplace=True )
		data_df_list = data_df.to_dict('records')
		try:
			self.col_insert_data('sg_bicluster_detail', data_df_list)
		except Exception as e:
			self.bind_object.logger.info("导入sg_bicluster_detail:%s信息出错:%s" % (bicluster_id, e,))
		else:
			self.bind_object.logger.info("导入sg_bicluster_detail:%s成功" % (bicluster_id,))





