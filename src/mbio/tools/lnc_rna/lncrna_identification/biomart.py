#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/4/10 10:09
@file    : biomart.py
"""

import json
import os
import unittest
from collections import defaultdict

from biocluster.agent import Agent
from biocluster.tool import Tool


# src/mbio/tools/lnc_rna/lncrna_identification/known_lncrna_identify.py


class BiomartAgent(Agent):
    def __init__(self, parent):
        super(BiomartAgent, self).__init__(parent)
        options = [
            {'name': 'biomart', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'biomart_type', 'type': 'string', 'default': 'type1'},
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
        if not self.option('biomart').is_set:
            raise Exception('biomart is not existence')
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '3G'


class BiomartTool(Tool):

    def get_biomart(self):
        """
            {
                'gene': {
                    'ENSG00000275206': {'gene_description': des, 'gene_name': symbol}
                },
                'transcript': {
                    'ENST00000399306': 'ENSG00000073712'
                }
            }
        :return:
        """
        res_dict = defaultdict(dict)
        des_type = self.option('biomart_type')
        for line_items in self.option('biomart').csv_reader(remove_header=False):
            gene_id = line_items[0]
            tran_id = line_items[1]
            if des_type == "type3":
                symbol = ""
                des = line_items[3]
            elif des_type == "type2":
                symbol = line_items[2]
                des = line_items[5]
            elif des_type == "type1":
                symbol = line_items[2]
                des = line_items[7]
            else:
                raise Exception('type error: [type format: type1, type2, type3]')

            if gene_id not in res_dict['gene']:
                # 'gene_description': des, 'gene_name': symbol, 'transcript_id': tran_id
                res_dict['gene'][gene_id] = {'gene_description': des, 'gene_name': symbol}
            res_dict['transcript'][tran_id] = gene_id

        return res_dict

    def output_file(self, bio_dict):
        detail_txt = os.path.join(self.output_dir, 'biomart.xls')
        detail_json = os.path.join(self.output_dir, 'biomart.json')

        gene_dict = bio_dict['gene']
        trans_dict = bio_dict['transcript']
        line_demo = '{trans_id}\t{gene_id}\t{gene_name}\t{gene_description}\n'

        with open(detail_txt, 'w') as txt_handler:
            txt_handler.write('transcript_id\tgene_id\tgene_name\tgene_description\n')
            for t_id, g_id in trans_dict.items():
                txt_handler.write(line_demo.format(trans_id=t_id, gene_id=g_id, **gene_dict[g_id]))

        with open(detail_json, 'w') as json_handler:
            json.dump(bio_dict, json_handler, indent=4)

    def run(self):
        super(BiomartTool, self).run()
        biomart_dic = self.get_biomart()
        self.output_file(bio_dict=biomart_dic)
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
                "id": "biomart_" + str(random.randint(1, 10000)) + '_' + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "lnc_rna.lncrna_identification.biomart",
                "instant": False,
                "options": dict(
                    biomart='/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/'
                            'Homo_sapiens/Ensemble_release_89/biomart/Homo_sapiens.GRCh38.biomart_gene.txt',
                    biomart_type='type1',
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()


    unittest.main()