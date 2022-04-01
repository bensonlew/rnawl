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



class LncrnaClassWorkflow(Workflow):
    """
    blast比对流程
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(LncrnaClassWorkflow, self).__init__(wsheet_object)
        TARGET_TYPE = ('db', 'genome', 'file')
        options = [
            {"name": "task_id", "type": "string"},
            {"name": "task_name", "type": "string", "default": ""},
            {"name": "main_id", "type": "string", "default": ""},
            {'name': 'update_info', 'type': 'string'},
            {"name": "lnc_file", "type": "infile", "format": "lnc_rna.common"},
            {'name': "lnc_type", "type": 'string', 'default': "bed"},
            {"name": "m_file", "type": "infile", "format": "lnc_rna.common"},
            {'name': "m_type", "type": 'string', 'default': "bed"},
            {'name': 'length_bidirection', 'type': 'int', 'default': 1000}
        ]


        '''
            {"name": "m_gtf", "type": "infile", "format": "lnc_rna.gtf"},
            {'name': 'm_gtf_id', 'type': 'string', 'default': None},
            {"name": "m_bed", "type": "infile", "format": "lnc_rna.common"},
            {'name': 'm_bed_id', 'type': 'string', 'default': None},
            {"name": "lnc_gtf", "type": "infile", "format": "lnc_rna.gtf"},
            {'name': 'lnc_gtf_id', 'type': 'string', 'default': None},
            {"name": "lnc_bed", "type": "infile", "format": "lnc_rna.common"},
            {'name': 'lnc_bed_id', 'type': 'string', 'default': None},
        '''
        print options
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def check_options(self):
        return True

    def run(self):
        '''
        if self.option("lnc_gtf").is_set:
            self.gtf2bed("lnc_gtf", "lnc_bed")
        if self.option("m_gtf").is_set:
            self.gtf2bed("m_gtf", "m_bed")
        '''

        self.run_lncrna_class()
        super(LncrnaClassWorkflow, self).run()

    def gtf2bed(self, gtf_opt, bed_opt):
        self.option(gtf_opt).to_bed()
        bed_path = os.path.split(self.option(gtf_opt).prop['path'])[0]
        bed = os.path.join(bed_path, os.path.split(self.option("gtf").prop['path'])[1] + ".bed")
        self.option(bed_opt, bed)


    def run_lncrna_class(self):
        self.lncrna_class = self.add_tool('tool_lab.lncrna_class')

        opts = {
            "length_bidirection": self.option('length_bidirection')
        }
        if self.option("lnc_type") != "bed":
            opts.update({
                "lncrna_gtf": self.option('lnc_file').prop['path'],
            })
        else:
            opts.update({
                "lncrna_bed": self.option('lnc_file').prop['path'],
            })
        if self.option("m_type") != "bed":
            opts.update({
                "mrna_gtf": self.option('m_file').prop['path'],
            })
        else:
            opts.update({
                "mrna_bed": self.option('m_file').prop['path'],
            })

        self.lncrna_class.set_options(opts)
        self.lncrna_class.on("end", self.set_output)
        self.lncrna_class.run()



    def set_db(self):
        api = self.api.api("tool_lab.lncrna_class")
        lncrna_detail  = self.lncrna_class.output_dir + '/lncRNA_classifications.xls.stat.xls'
        api.add_lncrna_class_detail(ObjectId(self.option("main_id")), lncrna_detail)


    def set_output(self):
        self.set_db()
        for file in os.listdir(self.lncrna_class.output_dir):
            CopyFile().linkfile(os.path.join(self.lncrna_class.output_dir, file), os.path.join(self.output_dir, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [r'.*.fa', 'fa', '序列文件', 0],
        ])
        super(LncrnaClassWorkflow, self).end()
