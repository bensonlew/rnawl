# -*- coding: utf-8 -*-
# __author__ = 'xuxi'
# last_modify:20210911
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

class Rocplot(ApiBase):
	def __init__(self, bind_object):
		super(Rocplot, self).__init__(bind_object)
		self._project_type = 'tool_lab'

	def add_rocplot(self, s3_output, main_id=None, project_sn='tool_lab', task_id='tool_lab', params=None):
		if not main_id:
			name = "Rocplot"+'_'
			time_now = datetime.datetime.now()
			name += time_now.strftime("%Y%m%d_%H%M%S")
			main_info = dict(
				project_sn=project_sn,
				task_id=task_id,
				version="v1",
				name=name,
				created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
				desc='Rocplot',
				params= params if params else "null",
				status="start",
			)
			main_id = ToolLab().create_db_table('sg_rocplot', [main_info])
		else:
			main_id = ObjectId(main_id)

		try:
			update_dict = {"status":"end", "s3_pdf":s3_output}
			self.update_db_record('sg_rocplot', {"_id": ObjectId(main_id)}, update_dict)
		except Exception as e:
			self.bind_object.logger.info("导入ROC_pdf出错:%s" % e)
		else:
			self.bind_object.logger.info("导入ROC_pdf成功")





