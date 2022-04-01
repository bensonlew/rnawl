# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
#  20210409

from biocluster.api.database.base import Base, report_check
import datetime
from bson.objectid import ObjectId
import json
import re
from types import StringTypes
from biocluster.config import Config
from bson.son import SON
from mbio.api.database.whole_transcriptome.api_base import ApiBase


class Pmlst(ApiBase):

	def __init__(self, bind_object):
		super(Pmlst, self).__init__(bind_object)
		self._project_type = "tool_lab"

	def add_mlst(self, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
		if main_id is None:
			name = "Pmlst" + "_"
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
				desc='Pmlst',
				params=params,
				status="start")
			main_id = self.create_db_table('pmlst', [main_info])
		else:
			main_id = ObjectId(main_id)

			try:
				self.update_db_record('pmlst', main_id, )
			except Exception as e:
				self.bind_object.logger.error("导入pmlst数据出错:%s" % e)
			else:
				self.bind_object.logger.info("导入pmlst数据成功")
		return main_id

	def add_mlst_detail(self, main_id, file1, file2, sample):
		with open(file1, "r") as f:
			lines = f.readlines()
			list1 = [{"field": "sample_name", "filter": "false", "sort": "false", "title": "sample name",
					"type": "string"},
					{"field": "st_id", "filter": "false", "sort": "false", "title": "ST",
					"type": "string"}, ]
			funs = lines[0].strip().split("\t")
			for i in funs[2:]:
				list1.append({"field": i, "filter": "false", "sort": "false", "title": i, "type": "string"})
			table1_dict = {"column": list1, "condition": {}}
			data_list = []
			for line in lines[1:]:
				lin = line.strip().split("\t")
				insert_data = {
					"pmlst_id": ObjectId(main_id),
					"sample_name": lin[0],
					"st_id": lin[1],
				}
				for i in range(2,len(lin)):
					insert_data[funs[i]] = lin[i]
			data_list.append(insert_data)
			self.create_db_table('pmlst_st', data_list)
		with open (file2, "r") as g:
			lines = g.readlines()
			data_list1 =[]
			for line in lines[1:]:
				lin = line.strip().split("\t")
				insert_data = {
					"pmlst_id": ObjectId(main_id),
					"sample_name": sample,
					"locus": lin[0],
					"ali_len": lin[1],
					"all_len": lin[2],
					"allele": lin[3],
					"gaps": lin[4],
					"coverage": lin[5],
					"identity": lin[6],
				}
				data_list1.append(insert_data)
			self.create_db_table('pmlst_detail', data_list1)
		table2_dict = {
			"column": [
				{"field": "sample_name", "filter": "false", "sort": "false", "title": "sample name", "type": "string"},
				{"field": "locus", "filter": "false", "sort": "false", "title": "Locus", "type": "string"},
				{"field": "identity", "filter": "false", "sort": "false", "title": "Identity(%)", "type": "float"},
				{"field": "coverage", "filter": "false", "sort": "false", "title": "Coverage(%)", "type": "float"},
				{"field": "ali_len", "filter": "false", "sort": "false", "title": "Alignment Length(bp)", "type": "int"},
				{"field": "all_len", "filter": "false", "sort": "false", "title": "Allele Length(bp)", "type": "int"},
				{"field": "gaps", "filter": "false", "sort": "false", "title": "Gaps", "type": "int"},
				{"field": "allele", "filter": "false", "sort": "false", "title": "Allele", "type": "string"}],
			"condition": {}}
		table2_info = json.dumps(table2_dict, sort_keys=False, separators=(',', ':'))
		self.bind_object.logger.info(table1_dict)
		table1_info = json.dumps(table1_dict, sort_keys=False, separators=(',', ':'))
		try:
			self.update_db_record('pmlst', main_id, main_id=main_id, detail_table=table2_info, stat_table=table1_info, status="end")
		except Exception as e:
			self.bind_object.logger.error("导入pmlst_detail数据出错:%s" % e)
		else:
			self.bind_object.logger.info("导入pmlst_detail数据成功")