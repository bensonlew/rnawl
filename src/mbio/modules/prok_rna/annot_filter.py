# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modifiy = modified 2018.05.26
# 注释结果过滤module

from biocluster.module import Module
import os
import shutil
import re
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
from mbio.packages.prok_rna.copy_file import CopyFile
import unittest

class AnnotFilterModule(Module):
    """
    根据参数筛选blast比对结果
    """
    def __init__(self, work_id):
        super(AnnotFilterModule, self).__init__(work_id)
        options = [
            {"name": "nr_evalue", "type": "float", "default": 1e-5},
            {"name": "nr_similarity", "type": "float", "default": 0},
            {"name": "nr_identity", "type": "float", "default": 0},
            {"name": "swissprot_evalue", "type": "float", "default": 1e-5},
            {"name": "swissprot_similarity", "type": "float", "default": 0},
            {"name": "swissprot_identity", "type": "float", "default": 0},
            {"name": "eggnog_evalue", "type": "float", "default": 1e-5},
            {"name": "eggnog_similarity", "type": "float", "default": 0},
            {"name": "eggnog_identity", "type": "float", "default": 0},
            {"name": "kegg_evalue", "type": "float", "default": 1e-5},
            {"name": "kegg_similarity", "type": "float", "default": 0},
            {"name": "kegg_identity", "type": "float", "default": 0},
            {"name": "pfam_evalue", "type": "float", "default": 1e-5},
            {"name": "blast_nr_xml", "type": "infile", "format": "prok_rna.blast_xml"},
            {"name": "blast_string_xml", "type": "infile", "format": "prok_rna.blast_xml"},
            {"name": "blast_eggnog_xml", "type": "infile", "format": "prok_rna.blast_xml"},
            {"name": "blast_kegg_xml", "type": "infile", "format": "prok_rna.blast_xml"},
            {"name": "blast_swissprot_xml", "type": "infile", "format": "prok_rna.blast_xml"},
            {"name": "pfam_domain", "type": "infile", "format": "prok_rna.common"},
            {"name": "db", "type": "string", "default": "nr,swissprot,kegg,eggnog,pfam,go"},
            {"name": "blast2go_annot", "type": "infile", "format": "prok_rna.blast2go_annot"}

        ]
        self.add_option(options)

        self.nr_filter = self.add_tool("prok_rna.annotation.filter_annot")
        self.swissprot_filter = self.add_tool("prok_rna.annotation.filter_annot")
        self.eggnog_filter = self.add_tool("prok_rna.annotation.filter_annot")
        self.kegg_filter = self.add_tool("prok_rna.annotation.filter_annot")
        self.pfam_filter = self.add_tool("prok_rna.annotation.filter_annot")
        self.go_filter = self.add_tool("prok_rna.annotation.filter_annot")
        self.filter_tools = []
        self.step.add_steps('nr_filter', 'swissprot_filter', 'kegg_filter', 'eggnog_filter', 'pfam_filter', 'go_filter')

    def check_options(self):
        for db in self.option("db").split(","):
            if db == "nr":
                if not self.option("blast_nr_xml").is_set:
                    raise OptionError("nr输入文件没有设置", code = "25000201")
                else:
                    self.filter_tools.append(self.nr_filter)
            if db == "swissprot":
                if not self.option("blast_swissprot_xml").is_set:
                    raise OptionError("swissprot输入文件没有设置", code = "25000202")
                else:
                    self.filter_tools.append(self.swissprot_filter)
            if db == "kegg":
                if not self.option("blast_kegg_xml").is_set:
                    raise OptionError("kegg输入文件没有设置", code = "25000203")
                else:
                    self.filter_tools.append(self.kegg_filter)
            if db == "eggnog":
                if not self.option("blast_eggnog_xml").is_set:
                    raise OptionError("eggnog输入文件没有设置", code = "25000204")
                else:
                    self.filter_tools.append(self.eggnog_filter)
            if db == "pfam":
                if not self.option("pfam_domain").is_set:
                    raise OptionError("pfam_domain输入文件没有设置", code = "25000205")
                self.filter_tools.append(self.pfam_filter)
            if db == "go":
                if not self.option("blast2go_annot").is_set:
                    raise OptionError("go输入文件没有设置", code = "25000206")
                self.filter_tools.append(self.go_filter)


    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run_nr_blast_filter(self):
        options = {
            'xml': self.option("blast_nr_xml"),
            'types': "xml",
            'evalue': self.option('nr_evalue'),
            'identity': self.option('nr_identity'),
            'similarity': self.option('nr_similarity')
        }
        self.nr_filter.on('start', self.set_step, {'start': self.step.nr_filter})
        self.nr_filter.on('end', self.set_step, {'end': self.step.nr_filter})
        self.nr_filter.set_options(options)
        self.nr_filter.run()

    def run_go_filter(self):
        options = {
            'blast2go_annot': self.option("blast2go_annot"),
            'types': "go",
            'evalue': self.option('nr_evalue'),
            'identity': self.option('nr_identity'),
            'similarity': self.option('nr_similarity')
        }
        self.go_filter.on('start', self.set_step, {'start': self.step.go_filter})
        self.go_filter.on('end', self.set_step, {'end': self.step.go_filter})
        self.go_filter.set_options(options)
        self.go_filter.run()


    def run_swissprot_blast_filter(self):
        options = {
            'xml': self.option("blast_swissprot_xml"),
            'types': "xml",
            'evalue': self.option('swissprot_evalue'),
            'identity': self.option('swissprot_identity'),
            'similarity': self.option('swissprot_similarity')
        }
        self.swissprot_filter.on('start', self.set_step, {'start': self.step.swissprot_filter})
        self.swissprot_filter.on('end', self.set_step, {'end': self.step.swissprot_filter})
        self.swissprot_filter.set_options(options)
        self.swissprot_filter.run()

    def run_eggnog_blast_filter(self):
        options = {
            'xml': self.option("blast_eggnog_xml"),
            'types': "xml",
            'evalue': self.option('eggnog_evalue'),
            'identity': self.option('eggnog_identity'),
            'similarity': self.option('eggnog_similarity')
        }
        self.eggnog_filter.on('start', self.set_step, {'start': self.step.eggnog_filter})
        self.eggnog_filter.on('end', self.set_step, {'end': self.step.eggnog_filter})
        self.eggnog_filter.set_options(options)
        self.eggnog_filter.run()

    def run_kegg_blast_filter(self):
        options = {
            'xml': self.option("blast_kegg_xml"),
            'types': "xml",
            'evalue': self.option('kegg_evalue'),
            'identity': self.option('kegg_identity'),
            'similarity': self.option('kegg_similarity')
        }
        self.kegg_filter.on('start', self.set_step, {'start': self.step.kegg_filter})
        self.kegg_filter.on('end', self.set_step, {'end': self.step.kegg_filter})
        self.kegg_filter.set_options(options)
        self.kegg_filter.run()

    def run_pfam_filter(self):
        options = {
            'hmm': self.option("pfam_domain"),
            'types': "hmm",
            'evalue': self.option('pfam_evalue'),
        }
        self.pfam_filter.on('start', self.set_step, {'start': self.step.pfam_filter})
        self.pfam_filter.on('end', self.set_step, {'end': self.step.pfam_filter})
        self.pfam_filter.set_options(options)
        self.pfam_filter.run()

    def run(self):
        super(AnnotFilterModule, self).run()
        self.on_rely(self.filter_tools, self.set_output)

        self.run_nr_blast_filter()
        self.run_eggnog_blast_filter()
        self.run_kegg_blast_filter()
        self.run_swissprot_blast_filter()
        self.run_pfam_filter()
        self.run_go_filter()


    def set_output(self):
        for db in self.option("db").split(","):
            if db == "nr":
                CopyFile().linkdir(self.nr_filter.output_dir, self.output_dir + "/nr")
            if db == "swissprot":
                CopyFile().linkdir(self.swissprot_filter.output_dir, self.output_dir + "/swissprot")
            if db == "eggnog":
                CopyFile().linkdir(self.eggnog_filter.output_dir, self.output_dir + "/eggnog")
            if db == "kegg":
                CopyFile().linkdir(self.kegg_filter.output_dir, self.output_dir + "/kegg")
            if db == "pfam":
                CopyFile().linkdir(self.pfam_filter.output_dir, self.output_dir + "/pfam")
            if db == "go":
                CopyFile().linkdir(self.go_filter.output_dir, self.output_dir + "/go")
        self.end()

    def end(self):
        repaths = [
            [".", "", "blast输出目录"],
            ["blast.xml", "xml", "blast xml输出结果文件"],
            ["blast_table.xls", "xls", "blast xls输出结果文件"]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        super(AnnotFilterModule, self).end()
