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



class OrfPredictWorkflow(Workflow):
    """
    在目标DNA序列中搜寻开发阅读框，可以输出每个ORF所在的区域，并翻译成对应的蛋白序列
    此工具可以为新预测的DNA序列查找潜在的蛋白编码区
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(OrfPredictWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "input", "type": "infile", "format": "sequence.fasta"},  # FASTA序列文件
            {"name": "train_num", "type": "int", "default": 3000},
            {"name": "min_len", "type": "int", "default": 50},
            {"name": "genetic_code", "type": "string", "default": "universal"},
            {"name": "minimal_length", "type": "int", "default": 50},
            # NR（一级分类）['Animal','Plant','Fungi','Protist','All']
            {"name": "nr_database", "type": "string", "default": "All"},  # nr库类型
            # NR（二级分类）动物['Reptilia','Mammalia','Invertebrate','Fishes','Aves','Amphibia']植物['Algae','Spermatophyta','OtherPlant']
            {"name": "plant", "type": "string", "default": None},  # nr库类型
            {"name": "animal", "type": "string", "default": None},  # nr库类型
            {"name": "nr_sub_database", "type": "string", "default": None},  # nr库类型
            {"name": "single_best_only", "type": "string", "default": "yes"},
            {"name": "nr", "type": "string", "default": "yes"},
            {"name": "swissprot", "type": "string", "default": "yes"},
            {"name": "pfam", "type": "string", "default": "yes"},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        print options
        self.add_option(options)
        self.revise_infiles()
        self.cds_predict = self.add_module("tool_lab.annot_orfpfam")
        self.diamond = self.add_module("tool_lab.annot_mapdb")

        self.set_options(self._sheet.options())

    def check_options(self):
        if not self.option("input").is_set:
            raise OptionError('请输入FASTA文件')
        if self.option("animal"):
            self.option("nr_sub_database", self.option("animal"))
        if self.option("plant"):
            self.option("nr_sub_database", self.option("plant"))
        return True

    def run(self):
        self.run_diamond()
        super(OrfPredictWorkflow, self).run()

    def run_diamond(self):
        self.logger.info("开始运行blast注释")
        blast_opts = {
            'query': self.option("input").prop['path'],
            'method': 'diamond',
        }

        db2nrdb = {
            'Animal': "metazoa",
            'Plant': 'viridiplantae',
            'Fungi': 'fungi',
            'Protist': 'protist',
            'Archaea': 'archaea',
            'Bacteria': 'bacteria',
            'Viruses': 'viruses',
            'All': 'nr'
        }
        if self.option("nr_sub_database") in ["Spermatophyta", "Reptilia", "OtherPlant", "Mammalia", "Invertebrate",
                                              "Fishes", "Aves", "Amphibia", "Algae"]:
            blast_opts.update(
                {
                    'nr_db': self.option("nr_sub_database"),
                }
            )
        else:
            blast_opts.update(
                {
                    'nr_db': db2nrdb[self.option("nr_database")],
                }
            )
        self.diamond.set_options(blast_opts)
        self.diamond.on('end', self.run_cds_predict, 'diamond')
        self.diamond.run()

    def run_cds_predict(self):
        self.logger.info("开始运行cds预测")
        # t2g2u = glob.glob(self.assemble_filter.output_dir + "/*filter.gene_trans_map")[0]
        nrxml = os.path.join(self.diamond.output_dir, "nr", "blast.xml")
        swissxml = os.path.join(self.diamond.output_dir, "swissprot", "blast.xml")
        opts = {
            "fasta": self.option("input").prop['path'],
            "genetic_code": self.option("genetic_code"),
            "blast_nr_xml": nrxml,
            "p_length": self.option("minimal_length"),
            "single_best_only" : self.option("single_best_only"),
            "blast_swissprot_xml": swissxml,
            "isoform_unigene": None
            # "isoform_unigene" : os.path.join(self.assemble.output_dir,"Trinity.filter_t2g2u"),
        }
        self.cds_predict.set_options(opts)
        self.cds_predict.on("end", self.set_output, "cds_predict")
        self.cds_predict.run()

    def set_db(self):
        orf_predict_api = self.api.api("tool_lab.orf_predict")
        orf_predict_api.add_orf_detail(self.option("main_id"), os.path.join(self.cds_predict.output_dir, "all_predicted.xls"))
        table_dict = {
            "column": [
                {"field": "protein_id", "title": "protein_id", "filter": "false", "sort": "false", "type": "string"},
                {"field": "seq_id", "title": "seq_id", "filter": "false", "sort": "false", "type": "string"},
                {"field": "position", "title": "position", "filter": "false", "sort": "false", "type": "string"},
                {"field": "cds_len", "title": "cds_len", "filter": "false", "sort": "false", "type": "int"},
                {"field": "type", "title": "type", "filter": "false", "sort": "false", "type": "string"}
            ],
            "condition": {}
        }
        table_info = json.dumps(table_dict, sort_keys=True, separators=(',', ':'))

        orf_predict_api.update_db_record('orf_predict',
                                         query_dict={"main_id": ObjectId(self.option("main_id"))},
                                         update_dict={'status': 'end',
                                                                                             'table_data': table_info})


    def set_output(self):
        self.set_db()
        for file in os.listdir(self.cds_predict.output_dir):
            CopyFile().linkfile(os.path.join(self.cds_predict.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "ORF查找结果文件",0],
            [r'.*.fa', 'xls', 'ORF查找结果文件', 0],
        ])
        super(OrfPredictWorkflow, self).end()
