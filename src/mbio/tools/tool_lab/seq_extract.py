# -*- coding: utf-8 -*-
# __author__ = 'fwy'
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import re
import pandas as pd
import unittest


class SeqExtractAgent(Agent):
    """
    Used for seq extract from fasta .

    usage:
    input:
    fasta:
    >scaffold1
    AGTGTTCTCTTGTTACATTTGTTTCTGCTCAAACACCATTAAACAACCGTCCTTCAAACGCCCCATCC
    >scaffold2
    TCCTCTTCTTTTTTTTTCAAGCAGAAGACGGCATACGAGATTAGTCTTGGTGACTGGAGTTCAGACG
    >scaffold3
    TGTGCTCTTCCGATCTGAATCCTTTCAACTTGTTGGGTTGGTAGGCACACACGCTAGGCGGGAGAG
    list:
    scaffold3     139         1638       rrn16      -
    scaffold4     39093     40592     rrn16     +
    scaffold4     43033     45850     rrn23     +
    scaffold1     88942     91759     rrn23      -
    scaffold4     45963     46065     rrn4.5    +
    scaffold1     88727     88829     rrn4.5     -
    scaffold4     46339     46461     rrn5        +
    scaffold1     88331     88453     rrn5        -
    out:
    >rrn16 scaffold3 139 1638 -
    TCTCATGGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCATGCCTAACACATGCAAGTGTC
    >rrn16 scaffold4 39093 40592 +
    TCTCATGGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCATGCCTAACACATGCAAGTCGG
    >rrn23 scaffold4 43033 45850 +
    TTCAAACGACAAGAGGCTTGCGGTGGATACCTAGGCACCCAGAGACGAAGAAGGGCGTGGTAAAC
    >rrn23 scaffold1 88942 91759 -
    TTCAAACGACAAGAGGCTTGCGGTGGATACCTAGGCACCCAGAGACGAAGAAGGGCGTGGTAAAC
    >rrn4.5 scaffold4 45963 46065 +
    TAAGGTCACGGCAAGACGAGCCGTTTATCACTACGATAGGTGCTAAGTGGAGGTGCAGTAATGTATG
    >rrn4.5 scaffold1 88727 88829 -
    TAAGGTCACGGCAAGACGAGCCGTTTATCACTACGATAGGTGCTAAGTGGAGGTGCAGTAATGTATG
    >rrn5 scaffold4 46339 46461 +
    GATATTCTGGTGTCCCAGGCGTAGAGGAACCACACCGATCCATCTCGAACTTGGTGGTGAAACTCTA
    >rrn5 scaffold1 88331 88453 -
    GATATTCTGGTGTCCCAGGCGTAGAGGAACCACACCGATCCATCTCGAACTTGGTGGTGAAACTCT

    """
    def __init__(self, parent):
        super(SeqExtractAgent, self).__init__(parent)
        options = [
            {"name": "fasta_file", "type": "infile", "format": "ref_rna_v2.fasta"},  # ??????????????????
            {"name": "list", "type": "infile","format":"ref_rna_v2.common" },  #??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????+?????????/-?????????????????????????????? ??????????????????????????????????????????tab??????????????????????????????
        ]
        self.add_option(options)
        self.step.add_steps("seq_extract")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.seq_extract.start()
        self.step.update()

    def stepfinish(self):
        self.step.seq_extract.finish()
        self.step.update()

    def check_options(self):
        """
        ????????????????????????
        :return:
        """
        if not self.option('fasta_file').is_set:
            raise OptionError('???????????????????????????????????????')
        if not self.option('list').is_set:
            raise OptionError('?????????????????????')
        return True

    def set_resource(self):
        """
        ??????????????????????????????????????????????????? self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(SeqExtractAgent, self).end()


class SeqExtractTool(Tool):
    def __init__(self, config):
        super(SeqExtractTool, self).__init__(config)
        self._version = "v1.0.1"
        self.python_path = '/miniconda2/bin/python'
        self.perl =  'miniconda2/bin/'
        self.tool_path=self.config.PACKAGE_DIR+"/tool_lab/seq_extract.py"

    def run(self):
        """
        ??????
        :return:
        """
        super(SeqExtractTool, self).run()
        self.run_seq_extract()
        self.set_output()
        self.end()

    def run_seq_extract(self):
        self.logger.info("?????????fasta????????????????????????")
        extract_cmd = '{} {} -i {} -l {} -o {} --detail '\
                       .format (self.python_path, self.tool_path, self.option("fasta_file").prop["path"],
                          self.option("list").prop["path"],"extracted.fasta")
        command = self.add_command("seq_extract", extract_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("??????seq_extract??????")
        else:
            self.set_error("??????seq_extract????????????!")
            return False



    def set_output(self):
        """
        ????????????????????????output???????????????
        :return:
        """
        self.logger.info("??????????????????")
        try:
                result=os.path.join(self.work_dir, "extracted.fasta")
                link = os.path.join(self.output_dir, "extracted.fasta")
                if os.path.exists(link):
                    os.remove(link)
                os.link(result, link)
        except Exception as e:
            self.logger.info("????????????????????????{}".format(e))


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "seq_extract" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.seq_extract",
            "instant": True,
            "options": dict(
                # fasta_file=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                fasta_file ="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/seq_cut/input/example.fa" ,
                list="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/seq_cut/input/seqpos.list"
                # min_len=500,
                # max_len=10000
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
