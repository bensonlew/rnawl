# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'
# __last_modify__ = '2021.07.19'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from Bio import SeqIO


class ToolNglocAgent(Agent):
    """
    DemoAgent:
    version 1.0
    """

    def __init__(self, parent):
        super(ToolNglocAgent, self).__init__(parent)
        options = [
            {"name": "gene_faa", "type": "infile", "format": "sequence.fasta", "required": True},
            {"name": "pfam", "type": "string"},
            {"name": "sample", "type": "string"},
        ]
        self.add_option(options)

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = "5G"


class ToolNglocTool(Tool):
    def __init__(self, config):
        super(ToolNglocTool, self).__init__(config)
        self.python = '/program/Python/bin/python'

    def run_tiqu(self):
        """
        description
        :return:
        """
        list1 =[]
        with open(self.option("pfam"), "r") as f:
            lines =f.readlines()
            for line in lines[1:]:
                lin =line.strip().split("\t")
                list1.append(lin[0])
        list2 = list(set(list1))
        list3 = []
        for seqrecord in SeqIO.parse(self.option("gene_faa").prop['path'], "fasta"):
            if seqrecord.id in list2:
                list3.append(seqrecord)
        SeqIO.write(list3,self.output_dir + "/"+ self.option("sample") + ".pfam.fasta","fasta")


    def run(self):
        super(ToolNglocTool, self).run()
        self.run_tiqu()
        self.end()