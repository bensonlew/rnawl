# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modifies 20180427

"""gbk"""

import os
import re
import types
from biocluster.workflow import Workflow
from bson import ObjectId
import datetime
from biocluster.config import Config


class GbkWorkflow(Workflow):
    """
    交互gbk
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GbkWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "specimen_id", "type": "string"},
            {"name": "gbk_detail_id", "type": "string"},
            {"name": "gbk_id", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "gbk_file", "type": "infile", "format": "gene_structure.gbk"},  #
            {"name": "des", "type": "string"},  #
            {"name": "source", "type": "string"},  #
            {"name": "organism", "type": "string"},  #
            {"name": "author", "type": "string"},  #
            {"name": "title", "type": "string"},  #
            {"name": "journal", "type": "string"},  #
            {"name": "seq_type", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.gbk = self.add_tool("bacgenome.rename_gbk")
        self.file_path = self._sheet.output

    def run_gbk(self):
        if self.option('seq_type') in ['Scaffold','scaffold']:
            options = {
                "sample_name":self.option('specimen_id'),
                "gbk":self.option('gbk_file'),
                "des":self.option('des'),
                "source": self.option('source'),
                "organism": self.option('organism'),
                "author": self.option('author'),
                "title": self.option('title'),
                "journal": self.option('journal'),
            }
        else:
            options = {
                "sample_name": self.option('specimen_id'),
                "gbk": self.option('gbk_file'),
                "des": self.option('des'),
                "source": self.option('source'),
                "organism": self.option('organism'),
                "author": self.option('author'),
                "title": self.option('title'),
                "journal": self.option('journal'),
                "seq_type": self.option('seq_type')
            }
        self.gbk.set_options(options)
        self.gbk.on('end', self.set_db)
        self.gbk.run()

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.run_gbk()
        super(GbkWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        result = self.gbk.output_dir
        api_gbk = self.api.api('bacgenome.gbk')
        self.logger.info('dao biao ')
        self.logger.info(self.option('gbk_detail_id'))
        api_gbk.add_gbk_detail(path= result , gbk_detail_id=self.option('gbk_detail_id'),sample_name =self.option('specimen_id'),
                              task_id=self.option('task_id'), update=False, new_path =self.file_path)
        self.end()

    def end(self):
        repaths = [
            [".", "", "gbk文件目录"],
        ]
        regexps = [
            ['.*\.gbk$', '', 'gbk文件'],
            ['.*_Chromosome.*\.gbk$', '', '原色体gbk文件'],
            ['.*_Plasmid.*\.gbk$', '', '质粒gbk文件'],
        ]
        sdir = self.add_upload_dir(self.gbk.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(GbkWorkflow, self).end()