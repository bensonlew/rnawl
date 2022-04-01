# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modifies 20200601

import os,shutil
import re
import types
from biocluster.workflow import Workflow
from bson import ObjectId
import datetime
from biocluster.file import download, exists
from mbio.packages.bac_comp_genome.common_function import get_fasta, get_sample_from_tree,link_dir

class CompareTreeWorkflow(Workflow):
    """
    物种树和基因树的比较
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CompareTreeWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "gene_tree", "type": "infile", "format": "graph.newick_tree"},
            {"name": "species_tree", "type": "infile", "format": "graph.newick_tree"},
            {"name": "species_type", "type": 'string'},
            {"name": "sample_list", "type": 'string'},
            {"name": "gene_type", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.compare_tree = self.add_tool("bac_comp_genome.compare_tree")
        self.file_path = self._sheet.output

    def run_comapretree(self):
        self.compare_tree.set_options({
            "species_tree": self.option("species_tree"),
            "gene_tree": self.option("gene_tree"),
        })
        self.compare_tree.on("end", self.set_output)
        self.compare_tree.run()

    def run(self):
        self.run_comapretree()
        super(CompareTreeWorkflow, self).run()

    def set_output(self):
        link_dir(self.compare_tree.output_dir, self.output_dir)
        self.set_db()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        api_path = self.api.api("bac_comp_genome.common_api")
        main_id = self.option("main_id")
        api_path.add_main_detail2(self.output_dir + '/compare_result.xls', 'comp_tree_detail', main_id, "specimen_id", has_head=True, main_name='comp_id')
        self.end()

    def end(self):
        repaths = [
            [".", "", "", 0],
        ]
        regexps = [
            [r'.*\.xls$', 'xls', '', 0],
            [r'.*\.stat$', 'stat', '', 0]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(CompareTreeWorkflow, self).end()