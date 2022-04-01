#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/4/10 11:23
@file    : gtf_stat.py
"""
import json
import os
import unittest
from collections import defaultdict

from biocluster.agent import Agent
from biocluster.tool import Tool


class GtfStatAgent(Agent):
    def __init__(self, parent):
        super(GtfStatAgent, self).__init__(parent)
        options = [
            {'name': 'gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
        ]
        self.add_option(options)
        self.step.add_steps("biomart")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.biomart.start()
        self.step.update()

    def step_end(self):
        self.step.biomart.finish()
        self.step.update()

    def check_options(self):
        if not self.option('gtf').is_set:
            raise Exception('gtf is not set')
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '3G'


class GtfStatTool(Tool):

    def gtf_stat(self):
        res_dict = defaultdict(dict)
        with self.option('gtf') as gtf_handler:
            for line_splits, _ in gtf_handler:
                attr_dic = line_splits[8]
                trans_id = attr_dic['transcript_id']
                sub_dic = res_dict[trans_id]
                if len(sub_dic) == 0:
                    sub_dic['gene_id'] = attr_dic['gene_id']
                    sub_dic['strand'] = line_splits[6]
                    sub_dic['length'] = 0
                if line_splits[2] == 'exon':
                    exon_len = line_splits[4] - line_splits[3] + 1
                    sub_dic['length'] += exon_len

        return res_dict

    def output_file(self, data_dict):
        # detail_txt = os.path.join(self.output_dir, 'gtf_statistics.xls')
        detail_json = os.path.join(self.output_dir, 'gtf_statistics.json')

        line_demo = '{transcript_id}\t{gene_id}\t{strand}\t{length}\n'

        # with open(detail_txt, 'w') as txt_handler:
        #     txt_handler.write('transcript_id\tgene_id\tstrand\tlength\n')
        #     for t_id, sub_dic in data_dict.items():
        #         txt_handler.write(line_demo.format(transcript_id=t_id, **sub_dic))

        with open(detail_json, 'w') as json_handler:
            json.dump(data_dict, json_handler, indent=4)

    def run(self):
        super(GtfStatTool, self).run()
        stat_dic = self.gtf_stat()
        self.output_file(data_dict=stat_dic)
        self.end()


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """

        def test(self):
            import random
            from mbio.workflows.single import SingleWorkflow
            from biocluster.wsheet import Sheet

            data = {
                "id": "gtf_statistics_" + str(random.randint(1, 10000)) + '_' + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "lnc_rna.lncrna_identification.gtf_stat",
                "instant": False,
                "options": dict(
                    gtf='/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/assemble/method-cufflinks/output/NewTranscripts/new_genes.gtf',
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()


    unittest.main()
