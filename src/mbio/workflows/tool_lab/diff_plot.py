# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
import pandas as pd
import glob
from biocluster.config import Config
import time
from biocluster.file import getsize, exists
from biocluster.file import download
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError


class DiffPlotWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(DiffPlotWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "diff_file", "type": "infile", "format": "ref_rna_v2.common"},
            #  A csv file contains three columns, column one for gene ID (no duplicated allowed), column two for fold change and column three for pvalue.',
            {"name": "top", "type": "int", "default": 5},
            {"name": "fc", "type": "float", "default": 1},
            {"name": "color", "type": "int", "default": 1},
            # color scheme
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'},
            {'name': 'project_task_id', 'type': 'string'},
            {'name': 'relate_id', 'type': 'string'},
            {"name": 'source', 'type': 'string', 'default': 'tool_lab'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.diff_plot")
        self.api = self.api.api("tool_lab.api_base")
        self.set_options(self._sheet.options())

    def run(self):
        if self.option('source') == 'project':
            self.file_path = self.check_file_path()
            diff_file = self.download_s3_file(self.file_path, 'diff_plot.txt')
        else:
            diff_file = self.option('diff_file').path
        self.run_tool(diff_file)
        super(DiffPlotWorkflow, self).run()

    def check_options(self):
        if self.option('source') == 'tool_lab':
            if not self.option('diff_file').is_set:
                raise OptionError('差异基因文件必须输入')
        return True

    def run_tool(self, diff_file):
        opts = {
            'diff_file': diff_file,
            'top': self.option('top'),
            'fc': self.option('fc'),
            'method': self.option('color'),
        }
        self.tool.set_options(opts)
        self.tool.on('end', self.set_output)
        self.tool.run()

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
            print('sleep 10s')
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

    def set_output(self):
        for file in os.listdir(self.tool.output_dir):
            if os.path.exists(os.path.join(self.output_dir, file)):
                os.remove(os.path.join(self.output_dir, file))
            os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        query_dict = {"main_id": ObjectId(self.option("main_id"))}
        pdf = os.path.join(self._sheet.output, "diff_plot.pdf")
        update_dict = {"output": pdf, "status": "end"}
        self.api.update_db_record('diff_plot', query_dict, update_dict)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "差异基因表达排序图结果文件",0],
            [r'diff_plot.pdf', 'pdf', '差异基因表达排序PDF图', 0],
            [r'diff_plot.png', 'png', '差异基因表达排序PNG图', 0],
        ])
        super(DiffPlotWorkflow, self).end()