# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# __last_modify__ = '2021.07.19'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from Bio import SeqIO


class ToolPrimerAgent(Agent):
    """
    DemoAgent:
    version 1.0
    """

    def __init__(self, parent):
        super(ToolPrimerAgent, self).__init__(parent)
        options = [
            {"name": "genome_fa", "type": "infile", "format": "sequence.fasta", "required": True},
            {"name": "gene_gff", "type": "string"},
            {"name": "downstream", "type": "int"},
            {"name": "upstream", "type": "int"},
            {"name": "seq_id_new", "type": "string"},
            {"name": "seq_id", "type": "string"},
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        pass

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = "5G"


class ToolPrimerTool(Tool):
    def __init__(self, config):
        super(ToolPrimerTool, self).__init__(config)
        self.python = '/miniconda2/bin/python'

    def run_tiqu(self):
        """
        description
        :return:
        """
        dict ={}
        with open(self.option("gene_gff")) as f:
            lines =f.readlines()
            for line in lines[1:]:
                lin =line.strip().split("\t")
                if lin[0] == self.option("seq_id"):
                    location = lin[1].split("_ORF")[0]
                    dict[location] =[lin[2],lin[3],lin[4]]
        if os.path.exists(self.output_dir+"/"+self.option("seq_id_new")+".fasta"):
            os.remove(self.output_dir+"/"+self.option("seq_id_new")+".fasta")
        with open(self.output_dir+"/"+self.option("seq_id_new")+".fasta","w") as g:
            for seqrecord in SeqIO.parse(self.option("genome_fa").prop['path'], "fasta"):
                if seqrecord.id in dict.keys():
                    seq_len =len(seqrecord.seq)
                    if dict[seqrecord.id][2] == "+":
                        start = int(dict[seqrecord.id][0])-int(self.option("downstream"))
                        end = int(dict[seqrecord.id][1])+int(self.option("upstream"))
                        if start <=0:
                            seq = str(seqrecord.seq[0:end])
                            g.write(">{}\n{}\n".format(self.option("seq_id_new"), seq))
                        else:
                            seq = str(seqrecord.seq[start-1:end])
                            g.write(">{}\n{}\n".format(self.option("seq_id_new"), seq))
                        if end >= seq_len:
                            seq = str(seqrecord.seq[start - 1:seq_len])
                            g.write(">{}\n{}\n".format(self.option("seq_id_new"), seq))
                        else:
                            seq = str(seqrecord.seq[start - 1:end])
                            g.write(">{}\n{}\n".format(self.option("seq_id_new"), seq))

                    elif dict[seqrecord.id][2] == "-":
                        start = int(dict[seqrecord.id][1])-int(self.option("downstream"))
                        end = int(dict[seqrecord.id][0]) + int(self.option("upstream"))
                        if start <=0:
                            seq1 = str(seqrecord.seq[0:end])
                            seq = self.dna_reverse(self.dna_complement(seq1))
                            g.write(">{}\n{}\n".format(self.option("seq_id_new"), seq))
                        else:
                            seq1 = str(seqrecord.seq[start-1:end])
                            seq = self.dna_reverse(self.dna_complement(seq1))
                            g.write(">{}\n{}\n".format(self.option("seq_id_new"), seq))
                        if end >= seq_len:
                            seq1 = str(seqrecord.seq[start - 1:seq_len])
                            seq = self.dna_reverse(self.dna_complement(seq1))
                            g.write(">{}\n{}\n".format(self.option("seq_id_new"), seq))
                        else:
                            seq1 = str(seqrecord.seq[start - 1:end])
                            seq = self.dna_reverse(self.dna_complement(seq1))
                            g.write(">{}\n{}\n".format(self.option("seq_id_new"), seq))

    def dna_complement(self, sequence):
        sequence = sequence.upper()
        sequence = sequence.replace('A', 't')
        sequence = sequence.replace('T', 'a')
        sequence = sequence.replace('C', 'g')
        sequence = sequence.replace('G', 'c')
        return sequence.upper()

    def dna_reverse(self, sequence):
        sequence = sequence.upper()
        return sequence[::-1]

    def run(self):
        super(ToolPrimerTool, self).run()
        self.run_tiqu()
        self.end()