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
from biocluster.file import download, exists


class BlastWorkflow(Workflow):
    """
    基因比对查询
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BlastWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "specimens", "type": "string"},  #
            {"name": "sequence", "type": "string"},
            {"name": "method", "type": "string"},
            {"name": "evalue", "type": "string", "default": "1e-5"},
            {"name": "max_target_num", "type": "string", "default": "10"},
            {"name": "w_size", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "nr_id", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.geneblast = self.add_module("bacgenome.gene_blast")

    def run_blast(self):
        self.down_load_files()
        self.seq_to_file()
        evalue = str(min(float(self.option('evalue')),0.999999))
        options = {
            'sample_dir': self.sample_dir,
            'method': self.option('method'),
            "evalue": evalue,
            "max_target_num": self.option('max_target_num'),
            "w_size": self.option('w_size'),
            "search_sequence": self.sequence_file,
        }
        self.geneblast.set_options(options)
        self.geneblast.on('end', self.set_db)
        self.geneblast.run()

    def down_load_files(self):
        self.sample_dir = os.path.join(self.work_dir, "sample_dir")
        bacgenome = Bacgenome()
        bacgenome._config = Config()
        task_id = self.option("task_id")
        if not os.path.exists(self.sample_dir):
            os.mkdir(self.sample_dir)
        #if self.option("method") in ["blastn", "blastx", "tblastx"]:
        if self.option("method") in ["blastn", "tblastn"]:
            sequence_type = "fnn"
        else:  #["blastp","blastx"]
            sequence_type = "faa"
        samples_path = bacgenome.get_genefile_bysample(task_id, self.option("specimens"), type=sequence_type)
        for each in samples_path.keys():
            newfile = os.path.join(self.sample_dir, each)
            if not exists(samples_path[each]):
                self.logger.info("download file {} not exists".format(samples_path[each]))
            download(samples_path[each], newfile)

    def seq_to_file(self):
        seq = self.option("sequence")
        self.sequence_file = os.path.join(self.work_dir, "search_sequence")
        with open(self.sequence_file, "w") as file:
            seq = seq.replace('&gt;','>')
            if '>' not in seq:
                seq = '>query\n' + seq
            file.write(seq)

    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.run_blast()
        super(BlastWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        api_blast = self.api.api('bacgenome.gene_blast')
        self.linkdir(self.geneblast.output_dir, self.output_dir)
        all_files = os.listdir(self.geneblast.output_dir)
        for each in all_files:
            blast_detail_file = os.path.join(self.geneblast.output_dir, each)
            api_blast.add_blast_detail(self.option('main_id'), blast_detail_file, self.option("nr_id"))
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
            [".", "", "基因比对查询分析"],
        ]
        regexps = [
            [r'.*.xls', 'xls', '基因比对结果'],
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(BlastWorkflow, self).end()




