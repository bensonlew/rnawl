# -*- coding: utf-8 -*-
# __author__ = 'ysh'

import os
import re
import types
from biocluster.workflow import Workflow
from biocluster.config import Config
from mainapp.models.mongo.bacgenome import Bacgenome
from bson import ObjectId
import datetime
from biocluster.file import download,exists


class OrthomclWorkflow(Workflow):
    """
    同源基因分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(OrthomclWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "specimens", "type": "string"},  #
            {"name": "pv_cutoff", "type": "float", "default": "1e-5"},  # P-Value or E-Value Cutoff
            {"name": "pi_cutoff", "type": "int", "default": 0, "max": 100, "min": 0},  # Percent Identity Cutoff
            {"name": "inflation", "type": "float", "default": 1.5},  # Markov Inflation Index
            {"name": "pmatch_cutoff", "type": "int", "default": 0, "max": 100, "min": 0},  # Percent Match Cutoff
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "task_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.tool = self.add_tool("bacgenome.orthomcl")
        self.faa_dir = os.path.join(self.work_dir, "faa_dir")

    def run_orthomcl(self):
        self.down_load_files()
        options = {
            'fasta_dir': self.faa_dir,
            'pv_cutoff': self.option('pv_cutoff'),
            "pi_cutoff": self.option('pi_cutoff'),
            "inflation": self.option('inflation'),
            "pmatch_cutoff": self.option('pmatch_cutoff'),
        }
        self.tool.set_options(options)
        self.tool.on('end', self.set_db)
        self.tool.run()

    def down_load_files(self):
        bacgenome = Bacgenome()
        bacgenome._config = Config()
        task_id = self.option("task_id")
        if not os.path.exists(self.faa_dir):
            os.mkdir(self.faa_dir)
        samples_path = bacgenome.get_genefile_bysample(task_id, self.option("specimens"))
        for each in samples_path.keys():
            newfile = os.path.join(self.faa_dir, each + ".faa")
            if not exists(samples_path[each]):
                self.logger.info("download file {} not exists".format(samples_path[each]))
            download(samples_path[each], newfile)

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.run_orthomcl()
        super(OrthomclWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.linkdir(self.tool.output_dir, self.output_dir)
        venn_file = os.path.join(self.output_dir, "orthogenes_venn.xls")
        stat_file = os.path.join(self.output_dir, "orthomcl_stat.xls")
        up_path = os.path.join(self._sheet.output, "orthomcl_stat.xls")
        api_orthomcl = self.api.api('bacgenome.orthomcl')
        api_orthomcl.add_homology_detail(self.option('main_id'), self.option("specimens"), stat_file, up_path)
        api_orthomcl.add_homology_venn(self.option('main_id'), venn_file)
        self.end()

    def linkdir(self, dirpath, newdir):
        """
        link一个文件夹下的所有文件到本module的output目录
        """
        allfiles = os.listdir(dirpath)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        self.logger.info(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])

    def end(self):
        repaths = [
            [".", "", "同源基因分析"],
            ["orthogenes_venn.xls", "xls", "同源基因venn图结果"],
            ["orthomcl_stat.xls", "xls", "同源基因结果文件"],
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        super(OrthomclWorkflow, self).end()
