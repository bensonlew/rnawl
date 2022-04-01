# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

from bson.objectid import ObjectId
import datetime
import os
import json
import types
from mbio.api.database.whole_transcriptome.api_base import ApiBase

class GtdbTk(ApiBase):
    def __init__(self, bind_object):
        super(GtdbTk, self).__init__(bind_object)
        self._project_type = 'tool_lab'

    def add_gtdbtk(self,project_sn='tool_lab', main_id=None, task_id='tool_lab', params=None):
        if main_id is None:
            name = "GTDB-TK" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='GTDB-TK',
                params=params if params else "null",
                status="start",
            )
            main_id = self.create_db_table('sg_gtdb', [main_info])
        else:
            main_id = ObjectId(main_id)
        return main_id

    def add_gtdb_detail(self, file, main_id):
        table_dict = {"column": [
            {"field": "sample", "filter": "false", "sort": "false", "title": "Sample name", "type": "string"},
            {"field": "classification", "filter": "false", "sort": "false", "title": "classification",
             "type": "string"},
            {"field": "fastani_reference", "filter": "false", "sort": "false", "title": "fastani reference",
             "type": "float"},
            {"field": "reference_radius", "filter": "false", "sort": "false", "title": "fastani reference radius",
             "type": "float"},
            {"field": "fastani_ani", "filter": "false", "sort": "false", "title": "fastani ani", "type": "float"},
            {"field": "fastani_af", "filter": "false", "sort": "false", "title": "fastani af", "type": "float"},
            {"field": "closest_placement_ref", "filter": "false", "sort": "false",
             "title": "closest placement reference", "type": "float"},
            {"field": "closest_placement_radius", "filter": "false", "sort": "false",
             "title": "closest placement radius", "type": "float"},
            {"field": "closest_placement_ani", "filter": "false", "sort": "false", "title": "closest placement ani",
             "type": "float"},
            {"field": "classification_method", "filter": "false", "sort": "false", "title": "classification method",
             "type": "string"},
            {"field": "note", "filter": "false", "sort": "false", "title": "note", "type": "string"},
            {"field": "other_related_references", "filter": "false", "sort": "false",
             "title": "other related references", "type": "string"},
            {"field": "aa_percent", "filter": "false", "sort": "false", "title": "aa percent", "type": "float"},
            {"field": "red_value", "filter": "false", "sort": "false", "title": "red_value", "type": "string"},
            {"field": "warnings", "filter": "false", "sort": "false", "title": "warnings", "type": "string"}],
            "condition": {}}
        data_list = []
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                insert_data = {
                    "gtdb_id": ObjectId(main_id),
                    "sample": lin[0].split(".fasta")[0],
                    "classification": lin[1],
                    "fastani_reference": lin[2],
                    "reference_radius": lin[3] if lin[3] =="N/A" else float(lin[3]),
                    "fastani_ani": lin[5] if lin[5] =="N/A" else float(lin[5]),
                    "fastani_af": lin[6] if lin[6] =="N/A" else float(lin[6]),
                    "closest_placement_ref": lin[7],
                    "closest_placement_radius": lin[8] if lin[8] =="N/A" else float(lin[8]),
                    "closest_placement_ani": lin[11] if lin[11] =="N/A" else float(lin[11]),
                    "classification_method": lin[13],
                    "note": lin[14],
                    "other_related_references": lin[15],
                    "aa_percent": float(lin[16]),
                    "red_value": lin[18],
                    "warnings": lin[19],
                }
                data_list.append(insert_data)
        self.create_db_table('sg_gtdb_detail', data_list)
        table_info = json.dumps(table_dict, sort_keys=True, separators=(',', ':'))
        self.update_db_record('sg_gtdb', main_id, status="end", main_id=main_id, table_data=table_info)