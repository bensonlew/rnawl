# -*- coding: utf-8 -*-
# __author__ = 'konghualei, 20170421'
from mainapp.controllers.project.meta_controller import MetaController
from mainapp.libs.param_pack import group_detail_sort
from bson import ObjectId
import datetime
import shutil
import re,os
import time
from biocluster.file import getsize, exists
from biocluster.file import download
from biocluster.workflow import Workflow
from biocluster.config import Config
from bson import SON


class ToolGoCircWorkflow(Workflow):
    """
    基因集venn图
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ToolGoCircWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "update_info", "type": "string"},
            {'name': 'diff_exp', 'type': 'string'},
            {"name": "go_file", "type": "string"},
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
        self.set_db()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        try:
            main_id = ObjectId(self.option('main_id'))
        except:
            pass
        main_info = dict(
            task_id=self.option('task_id'),
            project_type='prok_rna',
            params=self.option('params'),
            status="end",
            main_id=main_id,
            relate_name=self.option('relate_name'),
            tool_type=self.option('tool_type'),
            relate_id=main_id,
            enrich_path=self.s3_path_enrich,
            diff_path=self.s3_path_diff
        )
        collection_name = 'tool_thurl'
        project_type = 'tool_lab'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        main_id_ = db[collection_name].insert_one(SON(main_info)).inserted_id
        conn = db[collection_name]
        task_id = conn.find_one({'_id': main_id_})['task_id']
        conn.update({'_id': main_id_, "task_id": task_id}, {"$set": {'main_id': main_id_}}, upsert=True)
        super(ToolGoCircWorkflow, self).end()

    def set_db(self):
        api = self.api.api("prok_rna.tool_lab_api")
        self.s3_path_enrich = '{}/{}'.format(self._sheet.output, 'go_file.txt')
        self.s3_path_diff = '{}/{}'.format(self._sheet.output, 'diff_exp.txt')
        api.add_go_circ(self.option('main_id'))
        self.set_output()

    def set_output(self):
        for each in ['go_file.txt', 'diff_exp.txt']:
            if os.path.exists(os.path.join(self.output_dir, each)):
                os.remove(os.path.join(self.output_dir, each))
            os.link(os.path.join(self.work_dir, each), os.path.join(self.output_dir, each))
        self.end()



