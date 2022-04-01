# -*- coding: utf-8 -*-
# __author__ = 'xuxi_20210729'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import os
import glob
from Bio import SeqIO
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
from biocluster.file import download
from bson.objectid import ObjectId
import time
from biocluster.file import getsize, exists


class CircularbarWorkflow(Workflow):
    """
    Circularbar base R package Circularbar
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CircularbarWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "infile", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "main_id", "type": "string"},
            {"name": 'source', 'type': 'string', 'default': 'tool_lab'},
            {'name': 'project_task_id', 'type': 'string'},
            {'name': 'relate_id', 'type': 'string'},
            {'name': 'update_info', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.Circularbar = self.add_tool("tool_lab.circularbar")
        self.set_options(self._sheet.options())

    def run(self):
        if self.option('source') == 'project':
            self.file_path = self.check_file_path()
            self.infile_path = self.download_s3_file(self.file_path, 'circularbar.txt')
        self.run_Circularbar()
        super(CircularbarWorkflow, self).run()

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
        elif os.path.exists(path):
            to_path = path
        elif exists(path):
            download(path, to_path)
        else:
            self.set_error('file can not find %s', variables=(path), code = '13700502')
        return to_path

    def check_options(self):
        if self.option('source') == 'tool_lab':
            if not self.option("infile").is_set:
                raise OptionError("必须设置输入表达量文件")
            return True

    def run_Circularbar(self):
        if self.option('source') == 'project':
            circularbar_infile = self.infile_path
        if self.option('source') == 'tool_lab':
            circularbar_infile = self.option('infile').prop['path']
        opts = {
            "infile":circularbar_infile
        }
        self.Circularbar.set_options(opts)
        self.Circularbar.on('end', self.set_db)
        self.Circularbar.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        if os.path.exists(os.path.join(self.Circularbar.output_dir, 'circular_bar.pdf')):
            pdf_s3_path = os.path.join(self._sheet.output,"circular_bar.pdf")
        else:
            pdf_s3_path = None
        Circularbar_api = self.api.api('tool_lab.circularbar')
        Circularbar_api.add_circularbar(
            s3_output=pdf_s3_path,
            main_id=self.option('main_id'),
        )
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.Circularbar.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "环形柱状图",0],
            ['./circular_bar.pdf', '', '环形柱状图', 0],
        ])
        super(CircularbarWorkflow, self).end()

