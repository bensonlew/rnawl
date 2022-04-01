# -*- coding: utf-8 -*-
# __author__ = "wangzhaoyue,shicaiping"

import glob
import os
import unittest

from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.ref_rna_v2.trans_step import class_code_count, tidy_code_num
from mbio.packages.ref_rna_v2.trans_step import gene_trans_exon, count_trans_or_exons
from mbio.packages.ref_rna_v2.trans_step import step_count


class StatisticsAgent(Agent):
    def __init__(self, parent):
        super(StatisticsAgent, self).__init__(parent)
        options = [
            {"name": "assemble_method", "type": "string", "default": "stringtie"},
            {"name": "all_transcripts", "type": "infile", "format": "ref_rna_v2.fasta"},
            {"name": "add_code_merged", "type": "infile", "format": "ref_rna_v2.gtf"},
            {"name": "new_transcripts", "type": "infile", "format": "ref_rna_v2.gtf"},
            {"name": "new_genes", "type": "infile", "format": "ref_rna_v2.gtf"},
            {"name": "old_transcripts", "type": "infile", "format": "ref_rna_v2.gtf"},
            {"name": "old_genes", "type": "infile", "format": "ref_rna_v2.gtf"},
        ]
        self.add_option(options)

    def check_options(self):
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "8G"

    def end(self):
        super(StatisticsAgent, self).end()


class StatisticsTool(Tool):
    def __init__(self, config):
        super(StatisticsTool, self).__init__(config)

    def run(self):
        super(StatisticsTool, self).run()
        self.cal_step_count()
        self.cal_class_code_count()
        self.cal_tidy_code_num()
        self.meddle_with_gtf()
        self.set_output()
        self.end()

    def cal_step_count(self):
        for path in glob.glob(os.path.join(self.output_dir, "trans_count_stat_*.txt")):
            os.remove(path)
        for step in (200, 300, 600, 1000):
            step_count(fasta_file=self.option("all_transcripts").path,
                       fasta_to_txt=os.path.join(self.work_dir, "transcripts.length.txt"), group_num=10, step=step,
                       stat_out=os.path.join(self.output_dir, "trans_count_stat_{}.txt".format(step)))

    def cal_class_code_count(self):
        class_code_count(gtf_file=self.option("add_code_merged").path,
                         code_num_trans=os.path.join(self.work_dir, "code_num.txt"))

    def cal_tidy_code_num(self):
        tidy_code_num(all_transcripts_fa=self.option("all_transcripts").path,
                      raw_code_num_txt=os.path.join(self.work_dir, "code_num.txt"),
                      new_code_num_txt=os.path.join(self.output_dir, "code_num.txt"))

    def meddle_with_gtf(self):
        input_files = list()
        for name in ("new_transcripts", "new_genes", "old_transcripts", "old_genes"):
            merged_gtf = self.option(name).path
            method = self.option("assemble_method").lower()
            gene_trans_file = os.path.join(self.work_dir, "{}.trans".format(os.path.basename(merged_gtf)))
            trans_exon_file = os.path.join(self.work_dir, "{}.exon".format(os.path.basename(merged_gtf)))
            gene_trans_exon(merged_gtf, method, gene_trans_file, trans_exon_file)
            input_files.extend((gene_trans_file, trans_exon_file))
        for input_file in input_files:
            for step in (1, 5, 10, 20):
                count_file = "{}_{}.test.txt".format(input_file, step)
                final_file = os.path.join(self.output_dir, "{}_{}.txt".format(os.path.basename(input_file), step))
                count_trans_or_exons(input_file, step, count_file, final_file)

    def set_output(self):
        pass


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "statistics_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "ref_rna_v2.assembly.statistics",
            "instant": False,
            "options": dict(
                all_transcripts="/mnt/ilustre/users/sanger-dev/workspace/20200325/Single_new_transcripts_7392_8411"
                                "/NewTranscripts/output/all_transcripts.fa",
                add_code_merged="/mnt/ilustre/users/sanger-dev/workspace/20200325/Single_new_transcripts_7392_8411"
                                "/NewTranscripts/output/add_code_merged.gtf",
                new_transcripts="/mnt/ilustre/users/sanger-dev/workspace/20200325/Single_new_transcripts_7392_8411"
                                "/NewTranscripts/output/new_transcripts.gtf",
                new_genes="/mnt/ilustre/users/sanger-dev/workspace/20200325/Single_new_transcripts_7392_8411"
                          "/NewTranscripts/output/new_genes.gtf",
                old_transcripts="/mnt/ilustre/users/sanger-dev/workspace/20200325/Single_new_transcripts_7392_8411"
                                "/NewTranscripts/output/old_transcripts.gtf",
                old_genes="/mnt/ilustre/users/sanger-dev/workspace/20200325/Single_new_transcripts_7392_8411"
                          "/NewTranscripts/output/old_genes.gtf"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)
