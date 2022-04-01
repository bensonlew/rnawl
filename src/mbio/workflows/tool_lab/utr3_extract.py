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

class Utr3ExtractWorkflow(Workflow):
    """
    blast比对流程
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(Utr3ExtractWorkflow, self).__init__(wsheet_object)
        TARGET_TYPE = ('db', 'genome', 'file')
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "task_name", "type": "string", "default": ""},
            {"name": "main_id", "type": "string", "default": ""},
            {'name': 'update_info', 'type': 'string'},
            {"name": "genome_file", "type": "infile", "format": "ref_rna_v2.common"},
            {'name': "gene_list", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "gtf_file", "type": "infile", "format": "ref_rna_v2.common"},
            {'name': "type", "type": 'string', 'default': "genome"}, # genome or biomart
            {'name': 'species_name', 'type': 'string', 'default': ''},
            {'name': 'species_class', 'type': 'string', 'default': ''},
        ]

        print options
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def check_options(self):
        return True

    def run(self):
        self.run_utr3_extract()
        super(Utr3ExtractWorkflow, self).run()

    def gtf2bed(self, gtf_opt, bed_opt):
        self.option(gtf_opt).to_bed()
        bed_path = os.path.split(self.option(gtf_opt).prop['path'])[0]
        bed = os.path.join(bed_path, os.path.split(self.option("gtf").prop['path'])[1] + ".bed")
        self.option(bed_opt, bed)

    def run_utr3_extract(self):
        self.extract_utr3 = self.add_tool('tool_lab.extract_utr3')

        if self.option('type') == "genome":
            opts = {
                "type": self.option('type'),
                "ref": self.option('genome_file'),
                "gtf": self.option('gtf_file').prop['path']
            }
        else:
            species_class = self.get_species_class(self.option('species_name'))
            opts = {
                "type": self.option('type'),
                "species_name": self.option('species_name'),
                "species_class": species_class
            }

        if self.option("gene_list").is_set:
            opts.update({
                "gene_list": self.option("gene_list").prop['path']
            })

        self.extract_utr3.set_options(opts)
        self.extract_utr3.on("end", self.set_output)
        self.extract_utr3.run()

    def get_species_class(self, species_name):
        self.genome_db = self.api.api("gene_db.genome_db")
        ens_results = self.genome_db.search_ensembl2(species_name)
        if ens_results.count() < 1:
            self.set_error("未能找到对应基因组")
        else:
            e_dict = ens_results.next()
            # print e_dict
            if len(e_dict.get("division").split('Ensembl')) >= 2:
                spe_class = e_dict.get("division").split('Ensembl')[1].lower()
            else:
                # if spe_class == 'ENSEMBL':
                spe_class = "asia"
            return spe_class


    def set_db(self):
        pass
        # api = self.api.api("tool_lab.extract_utr3")
        # lncrna_detail  = self.extract_utr3.output_dir + '/extract_utr3ifications.xls.stat.xls'
        # api.add_extract_utr3_detail(ObjectId(self.option("main_id")), lncrna_detail)


    def set_output(self):
        self.set_db()
        for file in os.listdir(self.extract_utr3.output_dir):
            CopyFile().linkfile(os.path.join(self.extract_utr3.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [r'*.fa', 'fa', 'utr3序列文件', 0],
        ])
        super(Utr3ExtractWorkflow, self).end()
