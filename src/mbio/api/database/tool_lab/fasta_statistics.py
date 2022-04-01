# -*- coding: utf-8 -*-
# __author__ = 'zzg'
from bson.objectid import ObjectId
import datetime
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase


class FastaStatistics(ApiBase):
	def __init__(self, bind_object):
		super(FastaStatistics, self).__init__(bind_object)
		self._project_type = 'tool_lab'

	def add_FastaStatistics(self, fastastatistics_path, project_sn='tool_lab', main_id=None, task_id='tool_lab', params=None):
		if main_id is None:
			name = "FastaStatistics"+'_'
			time_now = datetime.datetime.now()
			name += time_now.strftime("%Y%m%d_%H%M%S")
			# if type(params) == dict:
			#	 params = json.dumps(params, sort_keys=True, separators=(',', ':'))
			main_info = dict(
				project_sn=project_sn,
				task_id=task_id,
				version="v1",
				name=name,
				created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
				desc='fastastatistics',
				params=[],
				status="start",
			)
			main_id = self.create_db_table('sg_fasta_statistics', [main_info])
		else:
			main_id = ObjectId(main_id)

		if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
			main_id = ObjectId(main_id)
		# add detail info
		table_dict = {"column": [{"field":"minimum","filter":"false","sort":"false","title":"Minimum","type":"int"},
										{"field":"longest","filter":"false","sort":"false","title":"Longest","type":"int"},
										{"field":"average","filter":"false","sort":"false","title":"average","type":"int"},
										{"field":"gc","filter":"false","sort":"false","title":"GC","type":"float"},
										{"field":"n50","filter":"false","sort":"false","title":"N50","type":"int"},
										{"field":"n90","filter":"false","sort":"false","title":"N90","type":"int"},
										{"field":"gap","filter":"false","sort":"false","title":"gap","type":"float"},
										{"field":"total_base_number","filter":"false","sort":"false","title":"Total_base_number","type":"int"},
										{"field":"contig_num","filter":"false","sort":"false","title":"contig_num","type":"int"}],
						"condition": {}}
		#table_data = json.dumps(table_dict)
		table_info = json.dumps(table_dict, sort_keys=True, separators=(',', ':'))
		with open(fastastatistics_path, "r") as f:
			lines = f.readlines()
			fastastatistics = lines[1].strip().split("\t")
			insert_data = {
				"statistics_id": main_id,
				"minimum": fastastatistics[0],
				"longest": fastastatistics[1],
				"average": fastastatistics[2],
				"gc": fastastatistics[3],
				"n50": fastastatistics[4],
				"n90": fastastatistics[5],
				"gap": fastastatistics[6],
				"total_base_number": fastastatistics[7],
				"contig_num": fastastatistics[8]
			}
		self.create_db_table('sg_fasta_statistics_detail', [insert_data])
		self.update_db_record('sg_fasta_statistics', main_id, status="end", main_id=main_id, table_data=table_info)
		return main_id
