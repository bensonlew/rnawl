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


class ToolVcfDistributionWorkflow(Workflow):
    """
    基因集venn图
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ToolVcfDistributionWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "update_info", "type": "string"},
            {'name': 'vcf_file', 'type': 'string'},
            {'name': 'tool_type', 'type': 'string'},
            {'name': 'ref_path', 'type': 'string', 'default': ''},  # obtain assembly level
            {"name": "main_id", "type": "string"},
            {"name": "relate_name", "type": "string"},
            {'name': 'params', 'type': 'string'}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run(self):
        self.start_listener()
        self.fire("start")
        if self.option('ref_path'):
            self.new_ref_path = os.path.join(self.config.SOFTWARE_DIR, 'database', self.option('ref_path').split('/app/database/')[1])
        else:
            self.new_ref_path = ''
        self.set_db()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        try:
            main_id = ObjectId(self.option('main_id'))
        except:
            pass
        main_info = dict(
            task_id=self.option('task_id'),
            project_type='ref_rna_v2',
            params=self.option('params'),
            status="end",
            main_id=main_id,
            relate_name=self.option('relate_name'),
            tool_type=self.option('tool_type'),
            relate_id=main_id,
            file_path=self.option('vcf_file')
        )
        if self.chrom_file:
            main_info.update({'chrom_path': self.chrom_file})
        collection_name = 'tool_thurl'
        project_type = 'tool_lab'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        main_id_ = db[collection_name].insert_one(SON(main_info)).inserted_id
        conn = db[collection_name]
        task_id = conn.find_one({'_id': main_id_})['task_id']
        conn.update({'_id': main_id_, "task_id": task_id}, {"$set": {'main_id': main_id_}}, upsert=True)
        super(ToolVcfDistributionWorkflow, self).end()

    def set_db(self):
        api = self.api.api("ref_rna_v3.tool_lab_api")
        if self.new_ref_path and os.path.exists(self.new_ref_path):
            self.chrom_file = self.check_chromosome()
        else:
            self.chrom_file = ''
        api.add_vcf_distribution(self.option('main_id'))
        self.set_output()

    def check_chromosome(self):
        chrom_file = os.path.join(self.output_dir, 'chrom_list.txt')
        with open(self.new_ref_path, 'r') as r, open(chrom_file, 'w') as c:
            for line in r:
                chrom, chr_type = line.strip().split('\t')
                if chr_type.lower() == 'chromosome':
                    c.write(str(chrom) + '\n')
        if os.path.getsize(chrom_file) > 0:
            chrom_out = '{}/{}'.format(self._sheet.output, 'chrom_list.txt')
        else:
            chrom_out = ''
        return chrom_out

    def set_output(self):
        self.end()



