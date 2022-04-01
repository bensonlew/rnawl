# -*- coding: utf-8 -*-

from bson import ObjectId
import datetime
import shutil
import re,os
import time
from biocluster.workflow import Workflow
from biocluster.config import Config
from bson import SON


class ToolPathwayNetworkWorkflow(Workflow):
    """
    基因集venn图
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ToolPathwayNetworkWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "enrich_id", "type": "string"},
            {"name": "anno_type", "type": "string"},
            {"name": "pvalue_padjust", "type": "string"},
            {"name": "params", "type": "string"},
            {'name': "tool_type", 'type': 'string', 'default': 'pathway_network'},
            {"name": "main_id", "type": "string"},
            {"name": "relate_name", "type": "string"},
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
            project_type='dia',
            params=self.option('params'),
            status="end",
            main_id=main_id,
            relate_name=self.option('relate_name'),
            tool_type=self.option('tool_type'),
            relate_id=main_id,
            file_path=self.s3_path
        )
        collection_name = 'tool_thurl'
        project_type = 'tool_lab'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        main_id_ = db[collection_name].insert_one(SON(main_info)).inserted_id
        conn = db[collection_name]
        task_id = conn.find_one({'_id': main_id_})['task_id']
        conn.update({'_id': main_id_, "task_id": task_id}, {"$set": {'main_id': main_id_}}, upsert=True)
        super(ToolPathwayNetworkWorkflow, self).end()

    def set_db(self):
        diff_plot = self.api.api("dia.tool_lab_api")
        if self.option("anno_type").lower() == "kegg":
            self.s3_path = '{}/{}'.format(self._sheet.output, 'kegg_enrich_table.txt')
        else:
            self.s3_path = '{}/{}'.format(self._sheet.output, 'go_enrich_table.txt')
        diff_plot.add_pathway_network(self.option('main_id'))
        self.set_output()

    def set_output(self):
        if self.option("anno_type").lower() == "kegg":
            diff_plot_path = os.path.join(self.work_dir, 'kegg_enrich_table.txt')
            if os.path.exists(os.path.join(self.output_dir, 'kegg_enrich_table.txt')):
                os.remove(os.path.join(self.output_dir, 'kegg_enrich_table.txt'))
            os.link(diff_plot_path, os.path.join(self.output_dir, 'kegg_enrich_table.txt'))
        else:
            diff_plot_path = os.path.join(self.work_dir, 'go_enrich_table.txt')
            if os.path.exists(os.path.join(self.output_dir, 'go_enrich_table.txt')):
                os.remove(os.path.join(self.output_dir, 'go_enrich_table.txt'))
            os.link(diff_plot_path, os.path.join(self.output_dir, 'go_enrich_table.txt'))
        self.end()


