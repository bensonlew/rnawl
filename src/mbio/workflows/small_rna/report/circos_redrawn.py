# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
import re
import os
from bson.objectid import ObjectId
from biocluster.config import Config
import pandas as pd
import glob
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class CircosRedrawnWorkflow(Workflow):
    """
    重绘circos工作流
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CircosRedrawnWorkflow, self).__init__(wsheet_object)
        options = [
            {'type': 'infile', 'name': 'config', 'format': 'small_rna.ini'},
            {'format': 'small_rna.fasta', 'type': 'infile', 'name': 'fasta'},
            {'format': 'small_rna.fasta', 'type': 'infile', 'name': 'ref'},
            {'default': '', 'type': 'string', 'name': 'gtf'},
            {'format': 'small_rna.common', 'type': 'infile', 'name': 'arf'},
            {'default': 10, "type": "int", 'name': "chr_num"},
            {'format': 'small_rna.common_dir', 'type': 'infile', 'name': 'map_stat'},
            {'default': None, "type": 'string', 'name': 'main_id'},
            {"name": "task_id", "type": "string", "default": None},
            {"name": "submit_location", "type": "string", "default": None},
            {"name": "samples", "type": "string", "default": None},
            {"name": "task_type", "type": "string", "default": None},
            {"name": "group_id", "type": "string", "default": None},
            {"name": "group_dict", "type": "string", "default": None},

            {'name': 'update_info', 'type': 'string'},
        ]

        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool_mapper_filter = self.add_tool("small_rna.mapper_stat_filter")

    def run(self):
        self.tool_mapper_filter.on("end", self.set_db)
        self.get_run_log()
        self.run_mapper_filter()
        super(CircosRedrawnWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("small_rna", table="sg_circos", main_id=self.option('main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def run_mapper_filter(self):
        options = dict(
            config=self.option("config"),
            chr_num=self.option("chr_num"),
            map_stat=self.option("map_stat"),
            arf=self.option("arf"),
            gtf=self.option("gtf"),
            ref=self.option("ref"),
            samples=self.option("samples")
        )

        self.tool_mapper_filter.set_options(options)
        self.tool_mapper_filter.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.logger.info('start set_output at {}'.format(self.__class__.__name__))
        for source in glob.glob(os.path.join(self.tool_mapper_filter.output_dir, '*')):
            basename = os.path.basename(source)
            link_name = os.path.join(self.output_dir, basename)
            os.link(source, link_name)
            self.logger.info('succeed in linking {} to {}'.format(source, link_name))
        self.logger.info('finish set_output at {}'.format(self.__class__.__name__))
        self.end()

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if re.match(r'tsanger:', workflow_output):
            workflow_output = workflow_output.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        elif re.match(r'sanger:', workflow_output):
            workflow_output = workflow_output.replace('sanger:', '/mnt/ilustre/data/')
        elif re.match(r'^\w+://\S+/.+$', workflow_output):
            pass
        else:
            self.set_error("json output wrong")
        self.workflow_output = workflow_output
        return workflow_output

    def end(self):
        target_dir = self.get_workflow_output_dir()
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.output_dir)

        result_dir.add_relpath_rules([
            [".", "", "circos结果目录", 0],
            ["./*.png", "png", "circos图片", 0],
            ["./*.svg", "svg", "circos图片", 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        db = Config().get_mongo_client(mtype="small_rna")[Config().get_mongo_dbname("small_rna")]
        col1 = db["sg_circos"]
        col1.update({"_id": ObjectId(self.option("main_id"))}, {"$set": {"result_dir": target_dir, "status":"end"}}, upsert=True)
        if self.option("samples") != "":
            col1.update({"_id": ObjectId(self.option("main_id"))}, {"$set": {"sample_list": self.option("samples").split(",")}}, upsert=True)
        super(CircosRedrawnWorkflow, self).end()
