# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import group_detail_sort
from bson import ObjectId
import datetime
import shutil
import re,os
import time
from biocluster.workflow import Workflow
from biocluster.config import Config
from bson import SON
from collections import OrderedDict
import json


class ToolCircularHeatmapWorkflow(Workflow):
    """
    基因集venn图
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ToolCircularHeatmapWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "update_info", "type": "string"},
            {'name': 'group_dict', 'type': 'string'},
            {"name": "metab_table", "type": "string"},
            {"name": 'ms_type', 'type': 'string'},
            {'name': 'tool_type', 'type': 'string'},
            {"name": "main_id", "type": "string"},
            {"name": "relate_name", "type": "string"},
            {'name': 'params', 'type': 'string'}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run(self):
        self.start_listener()
        self.fire("start")
        self.get_group_table()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        try:
            main_id = ObjectId(self.option('main_id'))
        except:
            pass
        main_info = dict(
            task_id=self.option('task_id'),
            project_type='metabolome_' + self.option('ms_type'),
            params=self.option('params'),
            status="end",
            main_id=main_id,
            relate_name=self.option('relate_name'),
            tool_type=self.option('tool_type'),
            relate_id=main_id,
            metab_path=self.s3_metab_path,
            group_path=self.s3_group_path,
        )
        collection_name = 'tool_thurl'
        project_type = 'tool_lab'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        main_id_ = db[collection_name].insert_one(SON(main_info)).inserted_id
        conn = db[collection_name]
        task_id = conn.find_one({'_id': main_id_})['task_id']
        conn.update({'_id': main_id_, "task_id": task_id}, {"$set": {'main_id': main_id_}}, upsert=True)
        super(ToolCircularHeatmapWorkflow, self).end()

    def set_db(self):
        api = self.api.api("metabolome.tool_lab_api")
        self.s3_metab_path = '{}/{}'.format(self._sheet.output, 'metab_table.txt')
        self.s3_group_path = '{}/{}'.format(self._sheet.output, 'group_table.txt')
        api.add_circular_heatmap(self.option('main_id'))
        self.set_output()

    def set_output(self):
        metab = os.path.join(self.work_dir, 'metab_table.txt')
        if os.path.exists(os.path.join(self.output_dir, 'metab_table.txt')):
            os.remove(os.path.join(self.output_dir, 'metab_table.txt'))
        os.link(metab, os.path.join(self.output_dir, 'metab_table.txt'))
        self.end()

    def get_group_table(self):
        group_dict = json.loads(self.option('group_dict'), object_pairs_hook=OrderedDict)
        group_file = os.path.join(self.output_dir, 'group_table.txt')
        with open(group_file, 'w') as g:
            g.write('#sample\tgroup\n')
            for group, ss in group_dict.items():
                for s in ss:
                    g.write(s + '\t' + group + '\n')
        self.set_db()

