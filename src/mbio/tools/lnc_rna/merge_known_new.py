# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.files.sequence.fastq import FastqFile
from mbio.files.sequence.file_sample import FileSampleFile
import os
import re
import unittest

class MergeKnownNewAgent(Agent):
    """
    用于workflow开始之前对所输入的文件进行详细的内容检测
    """

    def __init__(self, parent):
        super(MergeKnownNewAgent, self).__init__(parent)
        options = [
            {"name": "all_known_gtf", "type": "infile", "format": "lnc_rna.gtf"},
            {"name": "all_known_fa", "type": "infile", "format": "lnc_rna.fasta"},
            {"name": "new_mrna_fa", "type": "infile", "format": "lnc_rna.fasta"},
            {"name": "new_lncrna_fa", "type": "infile", "format": "lnc_rna.fasta"},
            {"name": "new_mrna_gtf", "type": "infile", "format": "lnc_rna.gtf"},
            {"name": "new_lncrna_gtf", "type": "infile", "format": "lnc_rna.gtf"},
            {"name": "all_known_list", "type": "infile", "format": "lnc_rna.common"},
            {"name": "new_mrna_list", "type": "infile", "format": "lnc_rna.common"},
            {"name": "new_lncrna_list", "type": "infile", "format": "lnc_rna.common"},
            {"name": "merge_gtf", "type": "outfile", "format": "lnc_rna.gtf"}
            ,
        ]
        self.add_option(options)
        self.step.add_steps("file_check")
        self.on('start', self.start_file_check)
        self.on('end', self.end_file_check)

    def start_file_check(self):
        self.step.file_check.start()
        self.step.update()

    def end_file_check(self):
        self.step.file_check.finish()
        self.step.update()

    def check_option(self):
        for opt in ["all_known_gtf", "all_known_fa", "new_mrna_fa", "new_lncrna_fa", "new_lncrna_fa", "new_mrna_gtf", "new_lncrna_gtf", "all_known_list"]:
            if not self.option(opt).is_set:
                raise OptionError("必须输入{}文件参数".format(opt))


    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = "5G"


class MergeKnownNewTool(Tool):
    """
    检查输入文件的格式是否符合要求
    """
    def __init__(self, config):
        super(MergeKnownNewTool, self).__init__(config)
        self.samples = list()

    def lnc_file(self):
        '''
        处理lncRNA相关数据库文件
        '''
        self.option("new_mrna_gtf").merge_gtf(self.option("new_lncrna_gtf").prop['path'], self.output_dir + "/new_all.gtf")
        self.option("new_mrna_fa").merge_fasta(self.option("new_lncrna_fa").prop['path'], self.output_dir + "/new_all.fa")

        self.option("all_known_gtf").merge_gtf(self.output_dir + "/new_all.gtf", self.output_dir + "/known_and_new.gtf")
        self.option("all_known_fa").merge_fasta(self.output_dir + "/new_all.fa", self.output_dir + "/known_and_new.fa")
        self.option("merge_gtf", self.output_dir + "/known_and_new.gtf")
        self.option("merge_gtf").to_isoform2unigene(self.output_dir + "/t2g2l.txt")
        os.system("cut -f 1,2 {} > {}".format(self.output_dir + "/t2g2l.txt", self.output_dir + "/t2g.txt"))


    def run(self):
        super(MergeKnownNewTool, self).run()
        self.lnc_file()
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "MergeKnownNew_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "lnc_rna.merge_known_new",
            "instant": False,
            "options": dict(
                all_known_gtf="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/lncrna/lncrna.gtf",
                all_known_fa="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/lncrna/lncrna.fa",
                new_mrna_gtf="/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_lnc_rna/lnc_predict_out/new_mrna.gtf",
                new_mrna_fa="/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_lnc_rna/lnc_predict_out/new_mrna.fa",
                new_lncrna_gtf="/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_lnc_rna/lnc_predict_out/novel_lncrna.gtf",
                new_lncrna_fa="/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_lnc_rna/lnc_predict_out/novel_lncrna.fa",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
