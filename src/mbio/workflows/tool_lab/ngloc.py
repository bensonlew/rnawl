# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import json
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError
from mbio.packages.lnc_rna.copy_file import CopyFile
from biocluster.file import getsize, exists
from biocluster.file import download
from biocluster.config import Config
import time


class NglocWorkflow(Workflow):
    """
    kegg 编辑
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(NglocWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "seq", "type": "infile", "format": "ref_rna_v2.common"},  # FASTA序列文件
            {"name": "go", "type": "infile", "format": "ref_rna_v2.common"},  # FASTA序列文件
            {'name': 'gram', 'type': 'string', 'default': "neg"},
            {'name': 'update_info', 'type': 'string'},
            {"name": 'source', 'type': 'string', 'default': 'tool_lab'},  # ['tool_lab', 'project'],
            {"name": 'relate_id', 'type': 'string'},
            {'name': 'task_id', 'type': 'string'},
            {'name': 'project_task_id', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        print options
        self.add_option(options)
        self.revise_infiles()
        self.multiloc = self.add_tool("tool_lab.multiloc")
        self.set_options(self._sheet.options())

    def check_options(self):
        if self.option('source') == 'tool_lab':
            if not self.option("seq").is_set:
                raise OptionError('请输入seq文件')
        return True

    def run(self):
        if self.option('source') == 'project':
            self.file_path = self.check_file_path()
            self.seq_path = self.download_s3_file(self.file_path, 'seq.fa')
        self.run_multiloc()
        super(NglocWorkflow, self).run()

    def run_multiloc(self):
        if self.option('source') == 'tool_lab':
            seq_file = self.option('seq').path
        elif self.option('source') == 'project':
            seq_file = self.seq_path
        self.logger.info("开始运行blast注释")
        opts = {
            "fa": seq_file,
            "gram": self.option("gram"),
            "species": "Bacteria"
        }

        # if self.option("go").is_set:
        #     opts.update({
        #         "go": self.option("go").prop[path]
        #     })
        self.multiloc.set_options(opts)
        self.multiloc.on('end', self.set_output)
        self.multiloc.run()


    def set_db(self):
        # pass
        multiloc_api = self.api.api("tool_lab.ngloc")
        multiloc_api.add_annotation_subloc_detail(self.option("main_id"), os.path.join(self.multiloc.output_dir, "multiloc.xls"))
        multiloc_api.add_annotation_subloc_bar(self.option("main_id"), os.path.join(self.multiloc.output_dir, "multiloc_stat.xls"))


        table_dict = {
            "column": [
                {"field": "accession_id", "title": "Accession", "filter": "false", "sort": "false", "type": "string"},
                {"field": "description", "title": "Description", "filter": "false", "sort": "false", "type": "string"},
                {"field": "subloc1", "title": "Subcellular Loc1", "filter": "false", "sort": "false", "type": "string"},
                {"field": "subloc1_prob", "title": "Loc1 Probability", "filter": "false", "sort": "false", "type": "string"},
                {"field": "subloc2", "title": "Subcellular Loc2", "filter": "false", "sort": "false", "type": "string"},
                {"field": "subloc2_prob", "title": "Loc2 Probability", "filter": "false", "sort": "false", "type": "string"},
                {"field": "subloc3", "title": "Subcellular Loc3", "filter": "false", "sort": "false", "type": "string"},
                {"field": "subloc3_prob", "title": "Loc3 Probability", "filter": "false", "sort": "false", "type": "string"}
            ],
            "condition": {}
        }
        table_info = json.dumps(table_dict, sort_keys=True, separators=(',', ':'))

        column_dict = {
            "name": "name",
            "data": "value",
            "condition": {"type": "column"}
        }
        column_data = json.dumps(column_dict, sort_keys=True, separators=(',', ':'))


        multiloc_api.update_db_record('ngloc',
                                      query_dict={"main_id": ObjectId(self.option("main_id"))},
                                      update_dict={'status': 'end',
                                                   'table_data': table_info,
                                                   'column_data': column_data
                                      }
        )

    def set_output(self):
        self.set_db()
        for file in os.listdir(self.multiloc.output_dir):
            CopyFile().linkfile(os.path.join(self.multiloc.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def check_file_path(self):
        collection_name = 'tool_thurl'
        project_type = 'tool_lab'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        conn_upset = db[collection_name]
        status = 'start'
        count_time = 0
        while status == 'start':
            if count_time > 600:
                self.set_error('超过十分钟还没有结果文件生成，请检查是否生成文件时报错')
                break
            time.sleep(10)
            print 'sleep 10s'
            try:
                upset = conn_upset.find_one(
                    {'task_id': self.option('project_task_id'), 'relate_id': ObjectId(self.option('relate_id'))})
                status = upset['status']
            except:
                pass
            count_time += 10
        upset = conn_upset.find_one(
            {'task_id': self.option('project_task_id'), 'relate_id': ObjectId(self.option('relate_id'))})
        file_path = upset['file_path']
        return file_path

    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        elif exists(path):
            download(path, to_path)
        else:
            self.set_error('file can not find %s', variables=(path), code = '13700502')
        return to_path

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "kegg图片编辑结果",0],
        ])
        super(NglocWorkflow, self).end()
