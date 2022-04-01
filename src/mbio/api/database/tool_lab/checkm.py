# -*- coding: utf-8 -*-
# __author__ = 'zzg'
from bson.objectid import ObjectId
import datetime
import json
import pickle
import unittest
import pandas as pd
from mbio.api.database.whole_transcriptome.api_base import ApiBase


class Checkm(ApiBase):
	def __init__(self, bind_object):
		super(Checkm, self).__init__(bind_object)
		self._project_type = 'tool_lab'

	def add_Checkm(self, Checkm_path, project_sn='tool_lab', main_id=None, task_id='tool_lab', params=None):
		if main_id is None:
			name = "Checkm"+'_'
			time_now = datetime.datetime.now()
			name += time_now.strftime("%Y%m%d_%H%M%S")
			main_info = dict(
				project_sn=project_sn,
				task_id=task_id,
				version="v1",
				name=name,
				created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
				desc='Checkm',
				params= params if params else "null",
				status="start",
			)
			main_id = self.create_db_table('sg_checkm', [main_info])
		else:
			main_id = ObjectId(main_id)

		if type(main_id) == str or type(main_id) == bytes or type(main_id) == unicode:
			main_id = ObjectId(main_id)
		table_dict = {"column": [{"field":"bin_id","filter":"false","sort":"false","title":"Bin_Id","type":"int"},
										{"field":"marker_lineage","filter":"false","sort":"false","title":"Marker_lineage","type":"int"},
										{"field":"genomes","filter":"false","sort":"false","title":"genomes","type":"int"},
										{"field":"markers","filter":"false","sort":"false","title":"markers","type":"float"},
										{"field":"marker_sets","filter":"false","sort":"false","title":"marker_sets","type":"int"},
										{"field":"a","filter":"false","sort":"false","title":"a","type":"int"},
										{"field":"b","filter":"false","sort":"false","title":"b","type":"int"},
										{"field":"c","filter":"false","sort":"false","title":"c","type":"int"},
										{"field":"d","filter":"false","sort":"false","title":"d","type":"int"},
										{"field":"e","filter":"false","sort":"false","title":"e","type":"int"},
										{"field":"f","filter":"false","sort":"false","title":"f","type":"int"},
										{"field":"completeness","filter":"false","sort":"false","title":"Completeness","type":"float"},
										{"field":"contamination","filter":"false","sort":"false","title":"Contamination","type":"int"},
										{"field":"strain_heterogeneity","filter":"false","sort":"false","title":"Strain_heterogeneity","type":"int"}],
						"condition": {}}
		table_info = json.dumps(table_dict, sort_keys=True, separators=(',', ':'))
		with open(Checkm_path, "r") as f:
			lines = f.readlines()
			Checkm = lines[1].strip().split("\t")
			insert_data = {
				"checkm_id": main_id,
				"bin_id": Checkm[0],
				"marker_lineage": Checkm[1],
				"genomes": Checkm[2],
				"markers": Checkm[3],
				"marker_sets": Checkm[4],
				"a": Checkm[5],
				"b": Checkm[6],
				"c": Checkm[7],
				"d": Checkm[8],
				"e": Checkm[9],
				"f": Checkm[10],
				"completeness": Checkm[11],
				"contamination": Checkm[12],
				"strain_heterogeneity": Checkm[13]
			}
		self.create_db_table('sg_checkm_detail', [insert_data])
		self.update_db_record('sg_checkm', main_id, status="end", main_id=main_id, table_data=table_info)
		return main_id
