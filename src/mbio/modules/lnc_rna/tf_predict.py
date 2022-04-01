#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import pandas as pd
import glob
import unittest
import shutil


class TfPredictModule(Module):
    """
    拆分后进行分别进行fasta预测
    """
    def __init__(self, work_id):
        super(TfPredictModule, self).__init__(work_id)
        options = [
            {'name': 's', 'type': 'string', 'default': 'plant', 'format': 'None'},
            {'name': 'organism', 'type': 'string', 'default': 'unknown', 'format': 'None'},
            {'name': 'blast_all', 'type': 'string', 'default': 'yes', 'format': 'None'},
            {'name': 'hmmscan', 'type': 'string',
             'default': '/mnt/ilustre/users/isanger/sg-users/litangjian/hmm/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan',
             'format': 'None'},
            {'name': 'hmmdb', 'type': 'string', 'default': '/mnt/ilustre/users/isanger/app/database/pfam_31/Pfam-A.hmm',
             'format': 'None'},
            {'name': 'seqfile', 'type': 'infile', 'default': '/mnt/ilustre/users/isanger/sg-users/deqing/TestFiles/TF/example.fa',
             'format': 'sequence.fasta'},
            {'name': 'E', 'type': 'float', 'default': '0.001', 'format': 'None'},
            {'name': 'domE', 'type': 'float', 'default': '0.0001', 'format': 'None'},
            {'name': 'cpu', 'type': 'int', 'default': '12', 'format': 'None'},
            {'name': 'tfdb', 'type': 'string', 'default': '/mnt/ilustre/users/isanger/app/database/TFDB/',
             'format': 'None'},
            {'name': 'diamond', 'type': 'string', 'default': '/mnt/ilustre/users/isanger/app/bioinfo/align/diamond-0.8.35/diamond',
             'format': 'None'},
            {'name': 'evalue', 'type': 'float', 'default': '0.0001', 'format': 'None'}, ]
        self.add_option(options)
        self.split_fasta = self.add_tool('lnc_rna.split_fasta')
        self.tools = list()

    def check_options(self):
        pass

    def run(self):
        super(TfPredictModule, self).run()
        self.run_split_fasta()

    def run_split_fasta(self):
        # split fasta
        split_options = dict(
            f=self.option('seqfile'),
            size=1000,
            prefix='split',
        )
        self.split_fasta.set_options(split_options)
        self.split_fasta.on('end', self.run_tf_predict)
        self.split_fasta.run()

    def run_tf_predict(self):
        # predict
        fasta_files = glob.glob(self.split_fasta.output_dir + '/split*fa')
        for fasta in fasta_files:
            tool_opts = dict(
                s=self.option('s'),
                organism=self.option('organism'),
                blast_all=self.option('blast_all'),
                seqfile=fasta,
                E=self.option('E'), # domE=self.option('domE'),
                evalue=self.option('evalue'),
            )
            tool = self.add_tool('lnc_rna.tf_predict')
            tool.set_options(tool_opts)
            self.tools.append(tool)
        if len(self.tools) == 1:
            self.tools[0].on("end", self.set_output)
            self.tools[0].run()
        else:
            self.on_rely(self.tools, self.set_output, "tf_predict")
        for tool in self.tools:
            tool.run()

    def set_output(self):
        tf_fa = list()
        blast_out = list()
        domain_out = list()
        final_tf_predict = list()
        for each_tool in self.tools:
            tmp = glob.glob(each_tool.output_dir + '/predicted_TFs.fa')
            if tmp:
                tf_fa.append(tmp[0])
            tmp = glob.glob(each_tool.output_dir + '/diamond.out.txt')
            if tmp:
                blast_out.append(tmp[0])
            tmp = glob.glob(each_tool.output_dir + '/domain_predict.txt')
            if tmp:
                domain_out.append(tmp[0])
            tmp = glob.glob(each_tool.output_dir + '/final_tf_predict.xls')
            if tmp:
                final_tf_predict.append(tmp[0])
        self.merge_result(tf_fa, self.output_dir+'/predicted_TFs.fa', has_header=False)
        self.merge_result(blast_out, self.output_dir+'/diamond.out.txt', has_header=False)
        self.merge_result(domain_out, self.output_dir + '/domain_predict.txt', has_header=True)
        self.merge_result(final_tf_predict, self.output_dir + '/final_tf_predict.xls', has_header=True)
        self.end()

    @staticmethod
    def merge_result(file_list, out_path, has_header=True):
        if not file_list:
            return
        with open(out_path, 'w') as fw:
            fr = open(file_list[0])
            fw.write(fr.read())
            fr.close()
            for each in file_list[1:]:
                fr = open(each)
                if has_header:
                    _ = fr.readline()
                fw.write(fr.read())
                fr.close()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([[".", "", "tf_predict"], ])
        super(TfPredictModule, self).end()

