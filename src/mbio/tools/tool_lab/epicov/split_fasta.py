# coding=utf-8
import os
import glob
import shutil
from Bio import SeqIO
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
__author__ = 'shicaiping'


class SplitFastaAgent(Agent):
    """
    split_fasta description
    """
    def __init__(self, parent):
        super(SplitFastaAgent, self).__init__(parent)
        options = [
            {'name': 'fasta', 'type': 'infile', 'default': '/mnt/ilustre/users/sanger-dev/sg-users/deqing/TestFiles/TF/all_pep.fa', 'format': 'ref_rna_v2.common'},
            {'name': 'num', 'type': 'int', 'default': 1}
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 2
        self._memory = "{}G".format('15')

    def end(self):
        super(SplitFastaAgent, self).end()


class SplitFastaTool(Tool):
    """
    split_fasta description
    """
    def __init__(self, config):
        super(SplitFastaTool, self).__init__(config)
        self.python_path = 'program/Python/bin/python'

    def run_split_fasta(self, fasta, prefix='split', num=1):
        i = 1
        j = 0
        seq_id = list()
        seq_sequence = list()
        for seq in SeqIO.parse(fasta, "fasta"):
            j += 1
            if j < num:
                seq_id.append(seq.id)
                seq_sequence.append(seq.seq)
            elif j == num:
                seq_id.append(seq.id)
                seq_sequence.append(seq.seq)
                subdir = os.path.join(self.output_dir, prefix + "_" + str(i))
                if os.path.exists(subdir):
                    shutil.rmtree(subdir)
                os.mkdir(subdir)
                with open(os.path.join(self.output_dir, subdir, prefix + "_" + str(i) + ".fa"), "w") as w:
                    for k, element in enumerate(seq_id):
                        w.write(">{}\n{}\n".format(seq_id[k], seq_sequence[k]))
                seq_id = list()
                seq_sequence = list()
                i += 1
                j = 0
        if seq_id:
            subdir = os.path.join(self.output_dir, prefix + "_" + str(i))
            if os.path.exists(subdir):
                shutil.rmtree(subdir)
            os.mkdir(subdir)
            with open(os.path.join(self.output_dir, subdir, prefix + "_" + str(i) + ".fa"), "w") as w:
                for k, element in enumerate(seq_id):
                    w.write(">{}\n{}\n".format(seq_id[k], seq_sequence[k]))

    def run(self):
        super(SplitFastaTool, self).run()
        self.run_split_fasta(self.option("fasta").prop['path'], num=self.option("num"))
        self.end()