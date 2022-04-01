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
from biocluster.config import Config



class TssSeqWorkflow(Workflow):
    """
    blast比对流程
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TssSeqWorkflow, self).__init__(wsheet_object)
        TARGET_TYPE = ('db', 'genome', 'file')
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "task_name", "type": "string", "default": ""},
            {"name": "is_db", "type": "string", "default": ""},
            {"name": "main_id", "type": "string", "default": ""},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'taxonmy', 'type': 'string', 'default': 'Animal'},
            {'name': 'ref_genome', 'type': 'string', 'default': None},
            {'name': 'genome_id', 'type': 'string', 'default': None},
            {'name': 'genome_id', 'type': 'string', 'default': None},

            {"name": "genome_file", "type": "infile", "format": "sequence.fasta"},  # FASTA序列文件
            {"name": "gtf_file", "type": "infile", "format": "ref_rna_v2.gtf"},
            {'name': 'genome_file_id', 'type': 'string', 'default': None},
            {'name': 'gtf_file_id', 'type': 'string', 'default': None},
            {"name": "id_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "id_file_id", "type": "string", 'default': None},
            {'name': 'up', 'type': 'int', 'default': 2000},
            {'name': 'down', 'type': 'int', 'default': 500}
        ]
        print options
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def check_options(self):
        return True

    def run(self):
        genome_path, gtf_path = self.get_genome_dir()
        self.run_tss_seq(genome_path, gtf_path)
        super(TssSeqWorkflow, self).run()

    def get_genome_dir(self):
        if self.option("genome_id"):
            database = Config().get_mongo_client(mtype='ref_rna_v2', dydb_forbid=True)[Config().get_mongo_dbname('ref_rna_v2', dydb_forbid=True)]
            collection = database['sg_genome_db']
            genome_info = collection.find_one({'genome_id': self.option('genome_id')})
            genome_path = os.path.join(self.config.SOFTWARE_DIR, 'database/Genome_DB_finish', genome_info["dna_fa"])
            gtf_path = os.path.join(self.config.SOFTWARE_DIR, 'database/Genome_DB_finish', genome_info["gtf"])
        elif self.option("genome_file").is_set and self.option("gtf_file").is_set:
            genome_path = self.option("genome_file").prop['path']
            gtf_path = self.option("gtf_file").prop['path']
        else:
            raise OptionError('基因序列未指定')
        return genome_path,gtf_path


    def run_tss_seq(self, ref, gtf):
        self.tss_seq = self.add_tool('tool_lab.tss_seq')

        opts = {
            "ref": ref,
            "gtf": gtf,
            "tss_up": self.option("up"),
            "tss_down": self.option("down"),
            'id_file': self.option('id_file')
        }
        self.tss_seq.set_options(opts)
        self.tss_seq.on("end", self.set_output)
        self.tss_seq.run()


    def set_db(self):
        sequence_tss_api = self.api.api("tool_lab.api_base")
        sequence_tss_api.update_db_record('tss_seq',
                                          query_dict={"main_id": ObjectId(self.option("main_id"))},
                                          update_dict={'status': 'end'})

    def set_output(self):
        self.set_db()
        for file in os.listdir(self.tss_seq.output_dir):
            CopyFile().linkfile(os.path.join(self.tss_seq.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [r'.*.fa', 'fa', '序列文件', 0],
        ])
        super(TssSeqWorkflow, self).end()
