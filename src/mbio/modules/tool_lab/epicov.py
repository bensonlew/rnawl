# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

from biocluster.module import Module
import sys
import os
import glob
import unittest
from Bio import SeqIO


class EpicovModule(Module):
    def __init__(self, work_id):
        super(EpicovModule, self).__init__(work_id)
        options = [
            {"name": "primer", "type": "infile", "format": "ref_rna_v2.common"},  # 引物序列文件
            {"name": "fasta_dir", "type": "infile", "format": "ref_rna_v2.common_dir"},  # 基因组序列文件
            {"name": "fasta", "type": "infile", "format": "ref_rna_v2.common"},  # 基因组序列文件
            {"name": "evalue", "type": "float", "default": 30000},
            {"name": "word_size", "type": "int", "default": 7},
            {"name": "max_hsps", "type": "int", "default": 10},  # Max targets per sequence
            {"name": "max_target_seqs", "type": "int", "default": 1000},  # Max targets to show
        ]
        self.add_option(options)
        self.split_fasta = self.add_tool("tool_lab.epicov.split_fasta")

    def check_options(self):
        self.logger.info('start check_options at {}'.format(self.__class__.__name__))
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run(self):
        super(EpicovModule, self).run()
        if self.option("fasta").is_set:
            self.logger.info("输入的单个fasta文件")
            self.run_split_fasta()
        elif self.option("fasta_dir").is_set:
            self.logger.info("输入的fasta文件夹")
            self.run_blast1()

    def run_split_fasta(self):
        options = {
            'fasta': self.option("fasta").prop['path'],
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

    def run_blast(self):
        fasta_files = glob.glob(self.split_fasta.output_dir + "/split_*/*.fa")
        self.blast_tools = list()
        if len(fasta_files) == 1:
            self.blast = self.add_tool("tool_lab.epicov.blast")
            self.blast_tools.append(self.blast)
            options = {
                'database': fasta_files[0],
                'query': self.option("primer").prop['path'],
                'evalue': self.option('evalue'),
                'word_size': self.option("word_size"),
                'max_hsps': self.option("max_hsps"),
                'max_target_seqs': self.option("max_target_seqs")
            }
            self.blast.set_options(options)
            self.blast.on("end", self.set_output)
            self.blast.run()
        else:
            for file in sorted(fasta_files):
                self.logger.info(file)
                self.blast = self.add_tool("tool_lab.epicov.blast")
                self.blast_tools.append(self.blast)
                options = {
                    'database': file,
                    'query': self.option("primer").prop['path'],
                    'evalue': self.option('evalue')
                }
                self.blast.set_options(options)
            self.on_rely(self.blast_tools, self.set_output)
            for tool in self.blast_tools:
                tool.run()

    def run_blast1(self):
        fasta_files = glob.glob(self.option("fasta_dir").prop['path'] + "/split_*/*.fa")
        self.blast_tools = list()
        if len(fasta_files) == 1:
            self.blast = self.add_tool("tool_lab.epicov.blast")
            self.blast_tools.append(self.blast)
            options = {
                'database': fasta_files[0],
                'query': self.option("primer").prop['path'],
                'evalue': self.option('evalue')
            }
            self.blast.set_options(options)
            self.blast.on("end", self.set_output)
            self.blast.run()
        else:
            for file in sorted(fasta_files):
                self.blast = self.add_tool("tool_lab.epicov.blast")
                self.blast_tools.append(self.blast)
                options = {
                    'database': file,
                    'query': self.option("primer").prop['path'],
                    'evalue': self.option('evalue')
                }
                self.blast.set_options(options)
            self.on_rely(self.blast_tools, self.set_output)
            for tool in self.blast_tools:
                tool.run()

    def set_output(self):
        primer_len = dict()
        for seq in SeqIO.parse(self.option("primer").prop['path'], "fasta"):
            primer_len[seq.id] = len(seq.seq)
        results = list()
        for tool in self.blast_tools:
            print(tool.output_dir)
            result = glob.glob(tool.output_dir + "/*.txt")[0]
            results.append(result)
        with open(os.path.join(self.output_dir, "total_result.txt"), "w") as w:
            for file in results:
                if os.path.getsize(file) > 0:
                    with open(file, "r") as f:
                        for line in f:
                            w.write(line)
        # default_header1 = ['Query-Name', 'Hit-Name', 'Identity-%', 'Length', 'Mismatch', 'Gapopen_num', 'Q-Begin',
        #                  'Q-End', 'Hsp-Begin', 'Hsp-End', 'Expect', 'Score', 'Hsp_qseq', 'Hsp_hseq', 'Match']
        with open(os.path.join(self.output_dir, "total_result.txt"), "r") as f:
            with open(os.path.join(self.output_dir, "target_result.txt"), "w") as w1, open(
                    os.path.join(self.output_dir, "target_result_align_detail.txt"), "w") as w2:
                target_list = list()
                flag = 0
                for line in f:
                    items = line.strip().split("\t")
                    if ' ' in items[14] and int(items[3]) == primer_len[items[0]]:
                        if items[1] in target_list:
                            continue
                        else:
                            target_list.append(items[1])
                        if flag == 0:
                            flag += 1
                            w1.write("\t".join(
                                ["引物探针", "毒株编号", "identity", "比对长度", "错配数", "gap数", "引物探针开始碱基位置", "引物探针结束碱基位置",
                                 "毒株序列开始碱基位置",
                                 "毒株序列结束碱基位置", "期望值", "比对得分"]) + "\n")
                            w1.write("\t".join(
                                [items[0], items[1], items[2], items[3], items[4], items[5], items[6], items[7],
                                 items[8],
                                 items[9], items[10], items[11]]) + "\n")
                            w2.write("\t".join([items[0], items[6], items[12], items[7]]) + "\n")
                            match = ""
                            for i, j in enumerate(items[14]):
                                if j == "|":
                                    match = match + "."
                                elif j == " ":
                                    match = match + items[13][i]
                            w2.write("\t".join([items[1], items[8], match, items[9]]) + "\n")
                        else:
                            w1.write("\t".join(
                                [items[0], items[1], items[2], items[3], items[4], items[5], items[6], items[7],
                                 items[8],
                                 items[9], items[10], items[11]]) + "\n")
                            match = ""
                            for i, j in enumerate(items[14]):
                                if j == "|":
                                    match = match + "."
                                elif j == " ":
                                    match = match + items[13][i]
                            w2.write("\t".join([items[1], items[8], match, items[9]]) + "\n")
                    else:
                        if items[1] not in target_list:
                            target_list.append(items[1])
                        else:
                            continue

        self.end()

    def end(self):
        super(EpicovModule, self).end()
