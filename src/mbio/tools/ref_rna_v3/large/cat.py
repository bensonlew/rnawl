# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
from Bio import SeqIO
import shutil
# import pandas as pd
__author__ = 'shicaiping'


class CatAgent(Agent):
    """
    split_fasta description
    """
    def __init__(self, parent):
        super(CatAgent, self).__init__(parent)
        options = [
            {'name': 'input_list', 'type': 'string', "default": ""},
            {'name': 'header', 'type': 'bool', "default": True},
            {'name': 'out_file', 'type': 'string', "default": ""},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        super(CatAgent, self).end()


class CatTool(Tool):
    """
    split_fasta description
    """
    def __init__(self, config):
        super(CatTool, self).__init__(config)

    def cat(self):
        input_list = self.option("input_list")
        output =self.option("out_file")
        flag = 0
        with open(input_list, "r") as f, open(output, "w") as w:
            for line in f:
                self.logger.info(line)
                flag += 1
                file = open(line.strip())
                if flag == 1:
                    for line1 in file:
                        w.write(line1)
                else:
                    if self.option("header"):
                        header = file.readline()
                    for line1 in file:
                        w.write(line1)


    def run(self):
        super(CatTool, self).run()
        self.cat()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "Cat" + str(random.randint(1, 10000)) + "_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v3.large.cat",
            "instant": False,
            "options": dict(
                input_list="/mnt/ilustre/users/isanger/workspace/20210315/Refrna_n34u_49qekvdhgok1ed33c1m5hf/CallSnpIndelSentieon/snp_annotation_list.txt",
                out_file="/mnt/ilustre/users/isanger/workspace/20210315/Refrna_n34u_49qekvdhgok1ed33c1m5hf/CallSnpIndelSentieon/snp_annotation_re.xls"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
