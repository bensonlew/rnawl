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



class SequenceBlastWorkflow(Workflow):
    """
    blast比对流程
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SequenceBlastWorkflow, self).__init__(wsheet_object)
        TARGET_TYPE = ('db', 'genome', 'file')
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "task_name", "type": "string", "default": ""},
            {"name": "query_file", "type": "infile", "format": "sequence.fasta"},  # FASTA序列文件
            {"name": "query_file_id", "type": "string"},
            {"name": "query_type", "type": "string", "default": "nucl"},
            {"name": "target_file", "type": "infile", "format": "sequence.fasta"},  # FASTA序列文件
            {"name": "target_file_id", "type": "string"},
            {"name": "target_type", "type": "string", "default": "nucl"},

            {"name": "is_target_file", "type": "string", "default": "false"},

            {"name": "actual_target_type", "type": "string", "default": "species"},
            {"name": "target_db", "type": "string", "default": None},

            {"name": "blast", "type": "string", "default": "blastx"},
            {"name": "num_alignment", "type": "int", "default": 5},
            {"name": "evalue", "type": "float", "default": 1e-5},
            {"name": "word_size", "type": "int", "default": 11},

            {'name': 'taxonmy', 'type': 'string', 'default': 'Animal'},
            {'name': 'species_name', 'type': 'string', 'default': None},
            {'name': 'genome_id', 'type': 'string', 'default': None},
            {'name': 'ref_db', 'type': 'string', 'default': None},
            {'name': 'db_id', 'type': 'string', 'default': None},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        print options
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self._sheet.options()

    def check_options(self):


        if self.option("query_type").lower() in ['prot', 'pro', 'protein']:
            self.option("query_type", 'prot')
        if self.option("query_type").lower() in ['nucl', 'nuc', 'dna']:
            self.option("query_type", 'nucl')
        if self.option("target_type").lower() in ['prot', 'pro', 'protein']:
            self.option("target_type", 'prot')
        if self.option("target_type").lower() in ['nucl', 'nuc', 'dna']:
            self.option("target_type", 'nucl')
        self.option('blast',  str(self.option('blast').lower()))
        blast_type = [
            ("nucl", "nucl", "blastn"),
            ("nucl", "nucl", "tblastx"),
            ("nucl", "prot", "blastx"),
            ("prot", "prot", "blastp"),
            ("prot", "nucl", "tblastn")
        ]

        self.map_blast = (self.option("query_type"), self.option("target_type"), str(self.option("blast").lower()))
        self.logger.info("map_blast {} is".format(self.map_blast))
        if self.map_blast in blast_type:
            pass
        else:
            raise OptionError('blast 软件选择错误')
        return True

    def run(self):
        db_path = self.get_database_dir()
        self.run_blast(db_path)
        super(SequenceBlastWorkflow, self).run()

    def get_database_dir(self):
        if self.option("is_target_file") == "true":
            db_path = self.option("target_file").prop['path']
        elif self.option("actual_target_type") == "species":
            database = Config().get_mongo_client(mtype='ref_rna_v2', dydb_forbid=True)[Config().get_mongo_dbname('ref_rna_v2', dydb_forbid=True)]
            collection = database['sg_genome_db']
            genome_info = collection.find_one({'genome_id': self.option('genome_id')})
            if self.option("target_type") == "nucl":
                db_path = os.path.join(self.config.SOFTWARE_DIR, 'database/Genome_DB_finish', genome_info["cds"])
            else:
                db_path = os.path.join(self.config.SOFTWARE_DIR, 'database/Genome_DB_finish', genome_info["pep"])

        elif self.option("actual_target_type") == "db":
            database = Config().get_mongo_client(mtype='ref_rna', ref=True)[Config().get_mongo_dbname('ref_rna', ref=True)]
            collection = database['blast_db']
            blastdb_info = collection.find_one({'name': self.option('target_db').lower()})
            db_path = Config().SOFTWARE_DIR + blastdb_info["path"]
        else:
            raise OptionError('数据库类型错误')
        return db_path


    def run_blast(self, db):
        self.blast_tool = self.add_tool('tool_lab.blast')

        opts = {
            "query_type": self.option("query_type"),
            "query": self.option("query_file"),
            "reference": db,
            "reference_type": self.option("target_type"),
            "outfmt": 6,
            "blast": self.option("blast"),
            "evalue": self.option("evalue"),

            "word_size": self.option("word_size"),
            "num_alignment": self.option("num_alignment"),

        }
        self.blast_tool.set_options(opts)
        self.blast_tool.on("end", self.set_output)
        self.blast_tool.run()


    def set_db(self):
        sequence_blast_api = self.api.api("tool_lab.sequence_blast")
        sequence_blast_api.add_blast_detail(self.option("main_id"), os.path.join(self.blast_tool.output_dir, "blast_result.xls"))
        show2field = {'Query-Name': 'query_name',
                      'Hit-Name': 'hit_name',
                      'Hit-Description': 'hit_description',
                      'HSP-Len': 'hsp-len',
                      'E-Value': 'e_value',
                      'Score': 'score',
                      'Identity-%': 'identity',
                      'Similarity-%': 'similarity'}
        table_dict = {
            "column": [
                {"field": v, "title": k, "filter": "false", "sort": "false", "type": "string"} for k,v in show2field.items()],
            "condition": {}
        }
        table_info = json.dumps(table_dict, sort_keys=True, separators=(',', ':'))

        sequence_blast_api.update_db_record('sequence_blast',
                                            query_dict={"main_id": ObjectId(self.option("main_id"))},
                                            update_dict={'status': 'end'})



        # 'table_data': table_info})


    def set_output(self):
        self.set_db()
        for file in os.listdir(self.blast_tool.output_dir):
            CopyFile().linkfile(os.path.join(self.blast_tool.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "ORF查找结果文件",0],
            [r'.*.xls', 'xls', 'ORF查找结果文件', 0],
        ])
        super(SequenceBlastWorkflow, self).end()
