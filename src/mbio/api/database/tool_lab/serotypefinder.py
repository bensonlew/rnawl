# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
#  20210430

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re
from types import StringTypes
from biocluster.config import Config
from bson.son import SON
from mbio.api.database.whole_transcriptome.api_base import ApiBase


class Serotypefinder(ApiBase):
	def __init__(self, bind_object):
		super(Serotypefinder, self).__init__(bind_object)
		self._project_type = "tool_lab"

	def add_mlst(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
		if main_id is None:
			name = "Serotypefinder" + "_"
			time_now = datetime.datetime.now()
			name += time_now.strftime("%Y%m%d_%H%M%S")
			if type(params) == dict:
				params = json.dumps(params, sort_keys=True, separators=(',', ':'))
			main_info = dict(
				project_sn=project_sn,
				task_id=task_id,
				version="v1",
				name=name,
				created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
				desc='Serotypefinder',
				params=params,
				status="start")
			main_id = self.create_db_table('serotypes', [main_info])
		else:
			main_id = ObjectId(main_id)

			try:
				self.update_db_record('serotypes', main_id, )
			except Exception as e:
				self.bind_object.logger.error("导入serotypes数据出错:%s" % e)
			else:
				self.bind_object.logger.info("导入serotypes数据成功")
		return main_id

	def add_serotypefinder_stat(self, main_id, file):
		data_list = []
		with open(file, "r") as f:
			lines = f.readlines()
			for line in lines[1:]:
				lin = line.strip().split("\t")
				insert_data = {
					"serot_id": ObjectId(main_id),
					"sample_name": lin[0],
					"o": lin[1],
					"o_gene": lin[2],
					"h": lin[3],
					"h_gene": lin[4],
					"s_type": lin[5],
				}
				data_list.append(insert_data)
		self.create_db_table('serotypes_st', data_list)
		table_dict = {
			"column": [
			{"field": "sample_name", "filter": "false", "sort": "false", "title": "sample name", "type": "string"},
			{"field": "o", "filter": "false", "sort": "false", "title": "O_type", "type": "string"},
			{"field": "o_gene","filter": "false","sort": "false","title": "O_type gene","type": "string"},
			{"field": "h","filter": "false","sort": "false","title": "H_type","type": "string"},
			{"field": "h_gene","filter": "false","sort": "false","title": "H_type gene","type": "string"},
			{"field": "s_type","filter": "false","sort": "false","title": "Serotype","type": "string"}],
			"condition":{}}
		table_info = json.dumps(table_dict, sort_keys=False, separators=(',', ':'))
		try:
			self.update_db_record('serotypes', main_id, main_id=main_id, stat_table=table_info, status="end")
		except Exception as e:
			self.bind_object.logger.error("导入serotypes_st数据出错:%s" % e)
		else:
			self.bind_object.logger.info("导入serotypes_st数据成功")


	def add_serotypefinder_detail(self, main_id, file, sample):
		data_list = []
		with open(file, "r") as f:
			lines = f.readlines()
			for line in lines[1:]:
				lin = line.strip().split("\t")
				insert_data = {
					"serot_id": ObjectId(main_id),
					"sample_name": sample,
					"location": lin[0],
					"serotype": lin[1],
					"gene": lin[2],
					"start": lin[3],
					"end": lin[4],
					"acc_num": lin[5],
					"identity": lin[6],
					"temp_len": lin[7],
					"hsp_len": lin[8],
				}
				data_list.append(insert_data)
		self.create_db_table('serotypes_detail', data_list)
		table_dict = {
			"column": [
				{"field": "sample_name", "filter": "false", "sort": "false", "title": "sample name", "type": "string"},
				{"field": "location", "filter": "false", "sort": "false", "title": "Location", "type": "string"},
				{"field": "serotype", "filter": "false", "sort": "false", "title": "Serotype", "type": "string"},
				{"field": "gene", "filter": "false", "sort": "false", "title": "Gene", "type": "string"},
				{"field": "start", "filter": "false", "sort": "false", "title": "Start", "type": "int"},
				{"field": "end", "filter": "false", "sort": "false", "title": "End", "type": "int"},
				{"field": "acc_num", "filter": "false", "sort": "false", "title": "AccNum", "type": "string"},
				{"field": "identity", "filter": "false", "sort": "false", "title": "Identity（%）", "type": "float"},
				{"field": "temp_len", "filter": "false", "sort": "false", "title": "Template Len（bp）", "type": "int"},
				{"field": "hsp_len", "filter": "false", "sort": "false", "title": "HSP Len（bp)", "type": "int"}],
			"condition": {}}
		table_info = json.dumps(table_dict, sort_keys=False, separators=(',', ':'))
		try:
			self.update_db_record('serotypes', main_id, main_id=main_id, detail_table=table_info, status="end")
		except Exception as e:
			self.bind_object.logger.error("导入serotypes_detail数据出错:%s" % e)
		else:
			self.bind_object.logger.info("导入serotypes_detail数据成功")