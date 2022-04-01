# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import os, glob, shutil
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from Bio import SeqIO
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError


class EpicovWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(EpicovWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "primer", "type": "infile", "format": "ref_rna_v2.common"},  # 引物序列文件
            {"name": "fasta", "type": "infile", "format": "ref_rna_v2.common"},  # 基因组序列文件
            {"name": "evalue", "type": "float", "default": 30000},
            {"name": "word_size", "type": "int", "default": 7},
            {"name": "max_hsps", "type": "int", "default": 10}, # Max targets per sequence
            {"name": "max_target_seqs", "type": "int", "default": 1000}, # Max targets to show
        ]
        self.add_option(options)
        self.revise_infiles()
        self.split_fasta = self.add_tool("tool_lab.epicov.split_fasta")
        self.split_primer = self.add_tool("tool_lab.epicov.split_fasta")
        self.set_options(self._sheet.options())

    def run(self):
        self.run_split_primer()
        self.run_split_fasta()
        super(EpicovWorkflow, self).run()

    def check_options(self):
        if not self.option('primer').is_set:
            raise OptionError('必须输入引物序列文件')
        if not self.option('fasta').is_set:
            raise OptionError('必须输入基因组序列文件')
        return True

    def run_split_fasta(self):
        num = self.option("max_target_seqs")/self.option("max_hsps")
        options = {
            'fasta': self.option("fasta").prop['path'],
            'num': num
        }
        self.split_fasta.set_options(options)
        self.split_fasta.on('end', self.build_index)
        self.split_fasta.run()

    def build_index(self):
        fasta_files = glob.glob(self.split_fasta.output_dir + "/split_*/*.fa")
        self.blast_index_tools = list()
        for file in sorted(fasta_files):
            self.build_index = self.add_tool("small_rna_v2.srna.blast_index")
            self.blast_index_tools.append(self.build_index)
            options = {
                'reference': file,
            }
            self.build_index.set_options(options)
        self.on_rely(self.blast_index_tools, self.run_blast)
        for tool in self.blast_index_tools:
            tool.run()

    def run_split_primer(self):
        options = {
            'fasta': self.option("primer").prop['path'],
            'num': 1
        }
        self.split_primer.set_options(options)
        self.split_primer.run()

    def run_blast(self):
        fasta_files = glob.glob(self.split_primer.output_dir + "/split_*/*.fa")
        if len(fasta_files) == 1:
            self.blast = self.add_module("tool_lab.epicov")
            self.blast_module.append(self.blast)
            options = {
                'primer': fasta_files[0],
                'fasta_dir': self.split_fasta.output_dir,
                'evalue': self.option('evalue'),
                'word_size': self.option("word_size"),
                'max_hsps': self.option("max_hsps"),
                'max_target_seqs': self.option("max_target_seqs")
            }
            self.blast.set_options(options)
            self.blast.on('end', self.set_output)
            self.blast.run()
        else:
            self.blast_module = list()
            for file in sorted(fasta_files):
                self.blast = self.add_module("tool_lab.epicov")
                self.blast_module.append(self.blast)
                options = {
                    'primer': file,
                    'fasta_dir': self.split_fasta.output_dir,
                    'evalue': self.option('evalue')
                }
                self.blast.set_options(options)
            self.on_rely(self.blast_module, self.set_output)
            for module in self.blast_module:
                module.run()

    def set_output(self):
        fasta_files = glob.glob(self.split_primer.output_dir + "/split_*/*.fa")
        for i, file in enumerate(sorted(fasta_files)):
            for seq in SeqIO.parse(file, "fasta"):
                seq_id = seq.id
                dir = os.path.join(self.output_dir, seq_id)
                if os.path.exists(dir):
                    shutil.rmtree(dir)
                os.mkdir(dir)
                for file in os.listdir(self.blast_module[i].output_dir):
                    os.link(os.path.join(self.blast_module[i].output_dir, file), os.path.join(dir, file))
        self.end()

    def end(self):
        super(EpicovWorkflow, self).end()
