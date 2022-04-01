# -*- coding: utf-8 -*-
# __author__ = 'xuxi'
# last_modify:20210908
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

class Mirnasearch(ApiBase):
	def __init__(self, bind_object):
		super(Mirnasearch, self).__init__(bind_object)
		self._project_type = 'tool_lab'

	def add_mirnasearch(self, result_file, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
		if not main_id:
			name = "Mirnasearch"+'_'
			time_now = datetime.datetime.now()
			name += time_now.strftime("%Y%m%d_%H%M%S")
			main_info = dict(
				project_sn=project_sn,
				task_id=task_id,
				version="v1",
				name=name,
				created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
				desc='Mirnasearch',
				params= params if params else "null",
				status="start",
			)
			main_id = ToolLab().create_db_table('sg_mirnasearch', [main_info])
		else:
			main_id = ObjectId(main_id)

		try:
			columns = ["target_entrez","database","mature_mirna_acc","target_ensembl","mature_mirna_id","target_symbol"]
			table_info = {'column': [{"filter": "false","field": i,"title": i.replace('_',' '),"type": 'string',"sort": "false"} for i in columns], 'condition': {}}
			update_dict = {"status":"end","table_data":json.dumps(table_info, sort_keys=True, separators=(',', ':'))}
			self.update_db_record('sg_mirnasearch', {"_id": ObjectId(main_id)}, update_dict)
			if os.path.exists(result_file):
				self.add_tsv_detail(result_file,mirnasearch_id=main_id)
		except Exception as e:
			self.bind_object.logger.info("导入miRNA搜索结果出错:%s" % e)
		else:
			self.bind_object.logger.info("导入miRNA搜索结果成功")

	def add_tsv_detail(self, tsv, mirnasearch_id):
		data_df = pd.read_table(tsv, header=0)
		data_df["mirnasearch_id"] = mirnasearch_id
		data_df_list = data_df.to_dict('records')
		try:
			self.col_insert_data('sg_mirnasearch_detail', data_df_list)
		except Exception as e:
			self.bind_object.logger.info("导入sg_mirnasearch_detail:%s信息出错:%s" % (mirnasearch_id, e,))
		else:
			self.bind_object.logger.info("导入sg_mirnasearch_detail:%s成功" % (mirnasearch_id,))




