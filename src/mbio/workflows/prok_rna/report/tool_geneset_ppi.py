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


class ToolGenesetPpiWorkflow(Workflow):
    """
    基因集venn图
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ToolGenesetPpiWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "submit_location", "type": "string"},
            {"name": "update_info", "type": "string"},
            {'name': 'gene_list', 'type': 'string'},
            {"name": "seq_dir", "type": "string"},
            {'name': 'tool_type', 'type': 'string'},
            {"name": "main_id", "type": "string"},
            {"name": "relate_name", "type": "string"},
            {'name': 'params', 'type': 'string'}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("prok_rna.tool_geneset_ppi")

    def run(self):
        # self.start_listener()
        # self.fire("start")
        # self.set_db()
        self.run_tool()
        super(ToolGenesetPpiWorkflow, self).run()

    def run_tool(self):
        if exists(os.path.join(self.option('seq_dir'), 'cds.fa')):
            seq_file = self.download_s3_file(os.path.join(self.option('seq_dir'), 'cds.fa'), 'cds.fa')
        elif exists(os.path.join(self.option('seq_dir'), 'Sequence_database', 'cds.fa')):
            seq_file = self.download_s3_file(os.path.join(self.option('seq_dir'), 'Sequence_database', 'cds.fa'), 'cds.fa')
        else:
            self.set_error("文件cds.fa不存在")
        options = {
            'geneset': self.option('gene_list'),
            'seq_file': seq_file,
        }
        self.tool.set_options(options)
        self.tool.on('end', self.set_db)
        self.tool.run()

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
            file_path=self.s3_path
        )
        collection_name = 'tool_thurl'
        project_type = 'tool_lab'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        main_id_ = db[collection_name].insert_one(SON(main_info)).inserted_id
        conn = db[collection_name]
        task_id = conn.find_one({'_id': main_id_})['task_id']
        conn.update({'_id': main_id_, "task_id": task_id}, {"$set": {'main_id': main_id_}}, upsert=True)
        super(ToolGenesetPpiWorkflow, self).end()

    def set_db(self):
        api = self.api.api("prok_rna.tool_lab_api")
        self.s3_path = '{}/{}'.format(self._sheet.output, 'geneset_file.txt')
        api.add_geneset_ppi(self.option('main_id'))
        self.set_output()

    def set_output(self):
        file_path = os.path.join(self.tool.work_dir, 'geneset_file.txt')
        if os.path.exists(os.path.join(self.output_dir, 'geneset_file.txt')):
            os.remove(os.path.join(self.output_dir, 'geneset_file.txt'))
        os.link(file_path, os.path.join(self.output_dir, 'geneset_file.txt'))
        self.end()

    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        if os.path.exists(path):
            to_path = path
        elif exists(path):
            download(path, to_path)
        else:
            self.set_error('file can not find %s', variables=(path,), code='13700502')
        return to_path


