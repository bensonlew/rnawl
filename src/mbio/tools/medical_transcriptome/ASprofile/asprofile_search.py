# coding=utf-8
import os
import unittest
import pandas as pd
import itertools
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.ref_rna_v2.functions import toolfuncdeco, runcmd
__author__ = 'zoujiaxun'


class AsprofileSearchAgent(Agent):
    """

    """
    def __init__(self, parent):
        super(AsprofileSearchAgent, self).__init__(parent)
        options = [

            dict(name="target_genes", type="string"),
            dict(name='as_result_merge', type='infile', format='ref_rna_v2.common'),
            dict(name='event_type', type='string'),
            dict(name='sample', type='string')

        ]
        self.add_option(options)

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        pass

    def set_resource(self):
        self._cpu = 1

        self._memory = '20G'

    def end(self):
        super(AsprofileSearchAgent, self).end()


class AsprofileSearchTool(Tool):
    def __init__(self, config):
        super(AsprofileSearchTool, self).__init__(config)

    def run(self):
        super(AsprofileSearchTool, self).run()
        self.run_asprofile_search()
        self.end()

    def run_asprofile_search(self):
        asprofile_detail = pd.read_table(self.option("as_result_merge").prop["path"])
        if not self.option("target_genes") == "all":
            target_genes =list()
            with open(self.option("target_genes")) as t:
                for i in t.readlines():
                    target_genes.append(i.strip())
        else:
            target_genes = "all"


        if  target_genes != "all":
            if self.option("event_type") == "" or self.option("event_type") == "all":
                search_df = asprofile_detail[(asprofile_detail["gene_id"].isin(target_genes)) & (
                                                     asprofile_detail["sample"] == self.option("sample"))]
                search_df.to_csv(os.path.join(self.output_dir, "target_detail.txt"), index=False, sep="\t")
            else:
                search_df = asprofile_detail[(asprofile_detail["event_type"] == self.option("event_type"))
                                             & (asprofile_detail["gene_id"].isin(target_genes)) & (
                                                         asprofile_detail["sample"] == self.option("sample"))]
                search_df.to_csv(os.path.join(self.output_dir, "target_detail.txt"), index=False, sep="\t")

        else:
            if self.option("event_type") == "" or self.option("event_type") == "all":
                search_df = asprofile_detail[(asprofile_detail["sample"] == self.option("sample"))]
                search_df.to_csv(os.path.join(self.output_dir, "target_detail.txt"), index=False, sep="\t")
            else:
                search_df = asprofile_detail[(asprofile_detail["event_type"] == self.option("event_type"))
                                             & (asprofile_detail["sample"] == self.option("sample"))]
                search_df.to_csv(os.path.join(self.output_dir, "target_detail.txt"), index=False, sep="\t")



    def set_output(self):
        # all ready write results to output
        pass




class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "ASprofile" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.ASprofile.asprofile",
            "instant": True,
            "options": {
                "transcripts": '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/ASprofile/test/merged.gtf',
                "hdrs": '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/ref_rna_v3/ASprofile/test/ref.fa.hdrs',
                'sample': 'Merge'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

