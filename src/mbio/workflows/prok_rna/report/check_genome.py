# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import os
from biocluster.config import Config
from bson import ObjectId
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.prok_rna.check_file import CheckFile


class CheckGenomeWorkflow(Workflow):
    """
    检查参考基因组准确性
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CheckGenomeWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "in_type", "type": "string"},
            {"name": "gtf", "type": "infile", "format": "gene_structure.gtf"},
            {"name": "genome", "type": "infile", "format": "prok_rna.common"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.check = self.add_tool('prok_rna.check_genome')

    def run_check(self):
        self.check.set_options({
            'in_type': self.option('in_type'),
            'in_file': self.option('gtf').prop['path'],
            'genome': self.option('genome').prop['path'],
        })
        self.check.on('end', self.end)
        self.check.run()

    def run(self):
        # self.start_listener()
        # self.fire("start")
        # if self.option("in_type").lower() == "gtf":
        #     check = CheckFile(self.option("genome").prop["path"], gtf=self.option("gtf").prop["path"],
        #                       work_dir=self.work_dir)
        # else:
        #     check = CheckFile(self.option("genome").prop["path"], gff=self.option("gtf").prop["path"],
        #                       work_dir=self.work_dir)
        # info = check.run()
        # self.logger.info(info)
        # if not isinstance(info, bool):
        #     self.set_error(info)
        self.run_check()
        super(CheckGenomeWorkflow, self).run()

    def end(self):
        super(CheckGenomeWorkflow, self).end()
