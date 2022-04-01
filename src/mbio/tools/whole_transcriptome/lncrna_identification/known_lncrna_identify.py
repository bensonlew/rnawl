#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/2/26 19:23
@file    : known_lncrna_identify.py
"""
import csv
import json
import os
import subprocess
import unittest
from collections import defaultdict
from multiprocessing import Process, Manager

import pandas as pd

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool

# src/mbio/tools/lnc_rna/lncrna_identification/known_lncrna_identify.py
from mbio.packages.lnc_rna.lnc_identification.predict_common_tools import CommandFactory


class KnownLncrnaIdentifyAgent(Agent):
    def __init__(self, parent):
        super(KnownLncrnaIdentifyAgent, self).__init__(parent)
        # lnc_db_gtf="/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/human/38/ensembl/lnc_db/Homo_sapiens.GRCh38.95.gtf",
        #                     lnc_db_fa="/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/human/38/ensembl/lnc_db/Homo_sapiens.GRCh38.95.gtf",
        #                     mrna_gtf='',
        #                     mrna_fa='',
        #                     exp_matrix="/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncrna_tools/transcript.tpm.matrix",
        options = [
            {'name': 'lnc_db_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
            {'name': 'lnc_db_fa', 'type': 'infile', 'format': 'lnc_rna.lnc_fasta'},
            {'name': 'ids_mapping', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'mrna_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
            {'name': 'mrna_fa', 'type': 'infile', 'format': 'lnc_rna.lnc_fasta'},
            {'name': 'exp_matrix', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'biomart_json', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'classify_info', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'database', 'type': 'string'}
        ]
        self.add_option(options)
        self.step.add_steps("known_lncrna_identify")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.known_lncrna_identify.start()
        self.step.update()

    def step_end(self):
        self.step.known_lncrna_identify.finish()
        self.step.update()

    def check_options(self):
        # super(KnownLncrnaIdentifyAgent, self).check_options()
        file_names = ['lnc_db_gtf', 'lnc_db_fa', 'mrna_gtf', 'mrna_fa', 'exp_matrix', 'ids_mapping']
        for name in file_names:
            if not self.option(name).is_set:
                raise Exception('%s is not existence' % name)
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        # self.add_upload_dir(self.output_dir)
        super(KnownLncrnaIdentifyAgent, self).end()


class KnownLncrnaIdentifyTool(Tool):

    def extr_classify(self):
        dic = {k: v for k, v in self.option('classify_info').csv_reader(remove_header=False)}
        return dic

    def extr_exp(self):
        ids_set = {i[0] for i in self.option('exp_matrix').csv_reader(remove_header=False)}
        return ids_set

    def extr_ids_dbs(self):
        self.logger.info("start read id table")
        dic = {d['transcript_id']: d['source'] for d in self.option('ids_mapping').csv_dict_reader()}
        self.logger.info("table dict is {}".format(dic))
        return dic

    def out_gtf(self, gtf_obj, out_file_path, needed_ids_set=None, out_ids_detail=False):
        all_ids = set()
        detail_dict = defaultdict(dict)
        with open(out_file_path, 'w') as out_handler, gtf_obj:
            for split_list, str_line in gtf_obj:
                attr_dict = split_list[8]
                trans_id = attr_dict.get('transcript_id')
                if trans_id is None:
                    continue
                all_ids.add(trans_id)

                if needed_ids_set and trans_id in needed_ids_set:
                    if out_ids_detail:
                        gene_id = attr_dict.get('gene_id')
                        if gene_id is None:
                            self.logger.debug('gene_id not in attribute')
                        feature = split_list[2]
                        if feature == 'exon':
                            sub_dic = detail_dict[trans_id]
                            if len(sub_dic) == 0:
                                sub_dic['strand'] = split_list[6]
                                sub_dic['lncrna_length'] = 0
                                sub_dic['gene_id'] = gene_id
                                sub_dic['lncrna_id'] = trans_id
                                sub_dic['chromosome'] = split_list[0]
                            sub_dic['lncrna_length'] += split_list[4] - split_list[3] + 1
                    out_handler.write(str_line)

        return all_ids, detail_dict

    def out_fa(self, fa_obj, out_file_path, needed_ids_set):
        # print needed_ids_set
        with open(out_file_path, 'w') as out_handler, fa_obj:
            for seq_id, seq_seq in fa_obj:
                #print "head is {}, str is {}".format(head, str_recored)
                #trans_id = head[0:].strip().split(' ')[0]
                #print trans_id
                if seq_id not in needed_ids_set:
                    continue
                out_handler.write(">" + seq_id + "\n" + str(seq_seq) + "\n")

    def out_detail(self, detail_dict, biomart_dict, classify_dic, mrna_detail_dict, lnc_dbs):
        """
        "_id" : ObjectId("5c99f0fe17b2bf12d2e521ef"),

    "database" : "ensembl",


        :param detail_dict:
        :param biomart_dict:
        :param classify_dic:
        :return:
        """
        trans_dic = biomart_dict['transcript']
        gene_dic = biomart_dict['gene']
        fields = ['lncrna_id', 'gene_id', 'chromosome', 'gene_name', 'lncrna_length', 'strand', 'gene_description',
                  'type', 'database']
        out_file = os.path.join(self.output_dir, 'known_lncrna_detail.xls')
        tmp_dic = {}
        database = self.option('database')
        all_lnc = set()
        with open(out_file, 'w') as out_handler:
            out_handler.write('\t'.join(fields) + '\n')
            line_demo = '\t'.join('{%s}' % i for i in fields) + '\n'
            for trans_id, dic in detail_dict.items():
                all_lnc.add(trans_id)
                # sub_gene_dic = gene_dic.get(trans_dic.get(trans_id), tmp_dic)
                sub_gene_dic = gene_dic.get(dic['gene_id'], tmp_dic)
                out_handler.write(line_demo.format(gene_name=sub_gene_dic.get('gene_name', ''),
                                                   gene_description=sub_gene_dic.get('gene_description', ''),
                                                   type=classify_dic[trans_id], database=lnc_dbs[trans_id], **dic))
                                                   # type=classify_dic.get(trans_id, '=='), database=lnc_dbs[trans_id], **dic))

        out_lnc_list = os.path.join(self.output_dir, 'known_lncrna_ids.list')
        with open(out_lnc_list, 'w') as out_handler:
            out_handler.write('\n'.join(all_lnc))

        out_lnc_list = os.path.join(self.output_dir, 'known_mrna_ids.list')
        with open(out_lnc_list, 'w') as out_handler:
            out_handler.write('\n'.join(i for i in mrna_detail_dict))

    def run(self):
        super(KnownLncrnaIdentifyTool, self).run()
        biomart = self.option('biomart_json').json_reader()
        classify_dic = self.extr_classify()
        extr_exp_ids = self.extr_exp()
        known_lncrna_gtf = os.path.join(self.output_dir, 'known_lncrna.gtf')
        all_known_lnc, lnc_detail_dic = self.out_gtf(self.option('lnc_db_gtf'), known_lncrna_gtf,
                                                     needed_ids_set=extr_exp_ids, out_ids_detail=True)
        known_lncrna_fa = os.path.join(self.output_dir, 'known_lncrna.fa')
        self.out_fa(self.option('lnc_db_fa'), known_lncrna_fa, needed_ids_set=extr_exp_ids)
        all_needed_mrna = extr_exp_ids - all_known_lnc
        known_mrna_gtf = os.path.join(self.output_dir, 'known_mrna.gtf')
        _, mrna_detail_dict = self.out_gtf(self.option('mrna_gtf'), known_mrna_gtf, needed_ids_set=all_needed_mrna,
                                           out_ids_detail=True)
        known_mrna_fa = os.path.join(self.output_dir, 'known_mrna.fa')
        self.out_fa(self.option('mrna_fa'), known_mrna_fa, needed_ids_set=all_needed_mrna)
        lnc_dbs = self.extr_ids_dbs()
        print "lnc_dbs is {}".format(lnc_dbs)
        self.out_detail(lnc_detail_dic, biomart, classify_dic, mrna_detail_dict, lnc_dbs)

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
                "id": "known_lncrna_identify_" + str(random.randint(1, 10000)) + '_' + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "lnc_rna.lncrna_identification.known_lncrna_identify",
                "instant": False,
                "options": dict(
                    lnc_db_gtf="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/lncrna/lncrna.gtf",
                    lnc_db_fa="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/lncrna/lncrna.fa",
                    ids_mapping="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/lncrna/ids_matrix.xls",
                    database='ensembl',
                    mrna_gtf='/mnt/ilustre/users/sanger-dev/workspace/20190402/LncDb_lnc_db_workflow_6536_6335/GtfFilter/output/mrna.gtf',
                    mrna_fa='/mnt/ilustre/users/sanger-dev/workspace/20190402/LncDb_lnc_db_workflow_6536_6335/GtfFilter/output/mrna.fa',
                    exp_matrix="/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncrna_tools/transcript.tpm.matrix",
                    biomart_json='/mnt/ilustre/users/sanger-dev/workspace/20190410/Single_biomart_2370_4297/Biomart/output/biomart.json',
                    classify_info='/mnt/ilustre/users/sanger-dev/workspace/20190411/Single_lncrna_classify_7366_5645/LncrnaClassify/output/known_lncRNA_classifications.xls'
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()


    unittest.main()
