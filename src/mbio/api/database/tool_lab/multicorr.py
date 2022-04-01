# -*- coding: utf-8 -*-
# __author__ = 'xuxi'
# last_modify:20210901
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

class Multicorr(ApiBase):
	def __init__(self, bind_object):
		super(Multicorr, self).__init__(bind_object)
		self._project_type = 'tool_lab'

	def add_multicorr(self, outfile_dir, s3_dir, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
		if not main_id:
			name = "Multicorr"+'_'
			time_now = datetime.datetime.now()
			name += time_now.strftime("%Y%m%d_%H%M%S")
			main_info = dict(
				project_sn=project_sn,
				task_id=task_id,
				version="v1",
				name=name,
				created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
				desc='Multicorr',
				params= params if params else "null",
				status="start",
			)
			main_id = ToolLab().create_db_table('sg_multicorr', [main_info])
		else:
			main_id = ObjectId(main_id)

		svg_list = [i[:-4] for i in os.listdir(outfile_dir) if ".svg" in i]
		svg_dict = {i[:-4]:os.path.join(s3_dir, i) for i in os.listdir(outfile_dir) if ".svg" in i}
		tsv_info = {'column': {}, 'condition': {}}
		if os.path.exists(os.path.join(outfile_dir ,"result_cor.txt")):
			data_df = pd.read_table(os.path.join(outfile_dir ,"result_cor.txt"), header=0)
			this_tsv_columns_ = data_df.columns.tolist()
			this_tsv_columns_a = [{"filter": "false","field": i,"title": i,"type": 'string',"sort": "false"} for i in this_tsv_columns_[:3]]
			this_tsv_columns_b = [{"filter": "false","field": i,"title": i,"type": 'float',"sort": "false"} for i in this_tsv_columns_[3:]]
			this_tsv_columns = this_tsv_columns_a + this_tsv_columns_b
			# this_tsv_columns.insert(0,{"filter": "false","field": 'gene',"title": 'gene',"type": 'string',"sort": "false"})
			tsv_info['column'] = this_tsv_columns
		try:
			update_dict = {
				"status":"end",
			 	"svg_list":json.dumps({"data":sorted(svg_list)}, sort_keys=True, separators=(',', ':')),
			 	"svg_dict":svg_dict,
			 	"detail_data":json.dumps(tsv_info, sort_keys=True, separators=(',', ':'))
			 }
			self.update_db_record('sg_multicorr', {"_id": ObjectId(main_id)}, update_dict)
		except Exception as e:
			self.bind_object.logger.info("导入组间相关性出错:%s" % e)
		else:
			self.bind_object.logger.info("导入组间相关性成功")

		if os.path.exists(os.path.join(outfile_dir ,"result_cor.txt")):
			self.add_tsv_detail(os.path.join(outfile_dir ,"result_cor.txt"),multicorr_id=main_id)

	def add_tsv_detail(self, tsv, multicorr_id):
		data_df = pd.read_table(tsv, header=0)
		data_df["multicorr_id"] = multicorr_id
		data_df_list = data_df.to_dict('records')
		try:
			self.col_insert_data('sg_multicorr_detail', data_df_list)
		except Exception as e:
			self.bind_object.logger.info("导入sg_multicorr_detail:%s信息出错:%s" % (multicorr_id, e,))
		else:
			self.bind_object.logger.info("导入sg_multicorr_detail:%s成功" % (multicorr_id,))

