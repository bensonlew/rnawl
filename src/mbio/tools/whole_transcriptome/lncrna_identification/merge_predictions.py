#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/2/28 23:00
@file    : merge_predictions.py
"""
import csv
import json
import os
import unittest
from itertools import chain
from pathlib import Path

from biocluster.agent import Agent
from biocluster.tool import Tool


class MergePredictionsAgent(Agent):
    def __init__(self, parent):
        super(MergePredictionsAgent, self).__init__(parent)
        # 'fasta_file', 'gtf_file', 'hexamer_dat', 'logit_model',
        # , 'transcript_len', 'cnci_score',
        # 'orf_len', 'cpat_score', 'taxonmy', 'exon_num'
        options = [
            {'name': 'predictions_dir', 'type': 'string'},
            {'name': 'new_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
            {'name': 'new_fasta', 'type': 'infile', 'format': 'lnc_rna.lnc_fasta'},
            {'name': 'biomart_json', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'classify_info', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'known_lnc_json', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'transcript2gene_info', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'pfam_mrna', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'tools', 'type': 'string'},
            {'name': 'identify_num', 'type': 'int', 'default': 1}
        ]
        self.add_option(options)
        self.step.add_steps("predictions_merge")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.predictions_merge.start()
        self.step.update()

    def step_end(self):
        self.step.predictions_merge.finish()
        self.step.update()

    def check_options(self):
        predictions_dir = Path(self.option('predictions_dir'))

        tools = self.option('tools').split(',')
        preds_files = [predictions_dir.joinpath(tool + '_output.txt') for tool in tools]
        for file_path in preds_files:
            if not file_path.is_file():
                raise Exception('%s is not existent' % file_path)
        if not self.option('new_gtf').is_set:
            raise Exception('%s is not existent' % self.option('new_gtf').prop['path'])
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(MergePredictionsAgent, self).end()


class MergePredictionsTool(Tool):
    """
    output_dir输出文件:
        new_lncrna_predict_detail.json: 新 lncRNA 预测详情表
        new_lncrna_stat.json: 新 lncRNA 预测统计详情表
        new_lncrna_ids.list: 新 lncRNA ids 列表
        new_lncrna.gtf: 已鉴定的新lncRNA GTF文件

        'biomart': os.path.join(self.biomart_tool.output_dir, 'biomart.json'),
        'classify_info': os.path.join(self.lnc_classify_tool.output_dir, 'all_new_rna_classifications.xls')
    """
    def classify(self):
        dic = {}
        for trans_id, class_type in self.option('classify_info').csv_reader(remove_header=False):
            dic[trans_id] = class_type
        return dic

    def file_reader(self, file, tool_name, final_dic, trans2gene_info, biomart, type_dict, extr_fields=None):
        """

        :param file:
        :param tool_name:
        :param final_dic:
        :param trans2gene_info: {'%transcript_id': {'strand': '', 'gene_id': '', 'length': 500}}
        :param biomart:
        :param extr_fields:
        :return:
        """
        temp_key = {'transcript_id', 'gene_id', 'label'}
        extr_fields = {k: tool_name + '_' + k for k in extr_fields if k not in temp_key}
        noncoding_list = []
        gene_info = biomart['gene']
        tmp_dict = {'gene_name': '', 'gene_description': ''}
        with open(file) as in_handler:
            for line_dic in csv.DictReader(in_handler, delimiter='\t'):
                t_id = line_dic['transcript_id']
                label = line_dic['label']
                label = 1 if label == 'noncoding' else 0

                sub_seq_detail = trans2gene_info.get(t_id)
                gene_id = sub_seq_detail['gene_id']
                sub_gene_info = gene_info.get(gene_id, tmp_dict)

                if t_id not in final_dic:
                    temp_dic = {'support_number': 0}
                    final_dic[t_id] = temp_dic
                    temp_dic['transcript_id'] = t_id
                    temp_dic['gene_name'] = sub_gene_info['gene_name']
                    temp_dic['gene_description'] = sub_gene_info['gene_description']
                    temp_dic['gene_id'] = gene_id
                    temp_dic['strand'] = sub_seq_detail['strand']
                    temp_dic['lncrna_length'] = sub_seq_detail['length']
                    temp_dic['type'] = type_dict.get(t_id, '')
                else:
                    temp_dic = final_dic[t_id]

                if label == 1:
                    noncoding_list.append((t_id, gene_id))
                temp_dic[tool_name + '_label'] = label
                temp_dic['support_number'] += label
                for name, new_name in extr_fields.items():
                    temp_dic[new_name] = line_dic[name]

        return list(zip(*noncoding_list))

    def trans2gene_info(self):
        with open(self.option('transcript2gene_info').path) as in_handler:
            return json.load(in_handler)

    def data_check(self, dic):
        tools = [i.strip() for i in self.option('tools').strip().split(',')]
        for tool in tools:
            t_score = tool + '_score'
            if t_score not in dic:
                dic[t_score] = 1111
                dic[tool + '_label'] = 0
        return dic

    def out_put(self, lncrna_detail_dic, lnc_stat_dic, tools):
        """
        "cpc_score": "0.104647",
        "pfam_score": 1111,
        "cpc_label": 1,
        "cnci_score": "-0.319488",
        "cpat_score": "0.00584430352671"

        :param lncrna_detail_dic:
        :param lnc_stat_dic:
        :return:
        """
        fields = ['transcript_id', 'gene_id', 'gene_name', 'strand', 'type']
        fields.extend((i + '_label' for i in tools))
        fields.extend((i + '_score' for i in tools if 'pfam' not in i))
        fields.extend(['gene_description', 'lncrna_length'])

        res_list = []
        identify_num = self.option('identify_num')
        new_lnc_list = set()
        detail_out_file = os.path.join(self.output_dir, 'novel_lncrna_predict_detail.xls')
        with open(detail_out_file, 'w') as out_handler:
            out_handler.write('\t'.join(fields) + '\n')
            line_demo = '\t'.join('{%s}' % i for i in fields) + '\n'
            for t_id, dic in lncrna_detail_dic.items():
                actual_num = dic.pop('support_number')
                if actual_num >= identify_num:
                    self.data_check(dic)
                    if 'pfam_score' in dic:
                        dic.pop('pfam_score')
                    new_lnc_list.add(t_id)
                    res_list.append(dic)
                    out_handler.write(line_demo.format(**dic))

        stat_out_file = os.path.join(self.output_dir, 'novel_lncrna_stat.json')
        with open(stat_out_file, 'w') as out_handler:
            json.dump(lnc_stat_dic, out_handler, indent=4)

        new_lnc_out_file = os.path.join(self.output_dir, 'novel_lncrna_ids.list')
        with open(new_lnc_out_file, 'w') as out_handler:
            out_handler.write('\n'.join(new_lnc_list))

        return new_lnc_list

    def get_new_mrna(self, pfam_predict):
        """
        'transcript_id\tlabel\n'
        :return:
        """
        new_mrnas = set()
        with open(pfam_predict) as in_handler:
            for dic in chain(self.option('pfam_mrna').csv_dict_reader(),
                             csv.DictReader(in_handler, delimiter='\t')):
                lable = dic['label']
                if lable == 'coding':
                    new_mrnas.add(dic['transcript_id'])
        return new_mrnas

    def out_new_gtf(self, new_lnc_list, known_lnc_set, new_mrna_set):
        novel_lnc_file = os.path.join(self.output_dir, 'novel_lncrna.gtf')
        novel_mrna_file = os.path.join(self.output_dir, 'novel_mrna.gtf')
        novel_mrna_ids = os.path.join(self.output_dir, 'novel_mrna_ids.list')
        new_lnc_set = {k for k in new_lnc_list}
        total_lncrans = set()
        novel_mrnas = set()
        with self.option('new_gtf') as gtf_handler, \
                open(novel_lnc_file, 'w') as lnc_handler, \
                open(novel_mrna_file, 'w') as mrna_handler, \
                open(novel_mrna_ids, 'w') as mrna_ids_handler:
            for split_line, str_line in gtf_handler:
                attr_dict = split_line[8]
                transcript_id = attr_dict.get('transcript_id', None)
                total_lncrans.add(transcript_id)
                if transcript_id is None or transcript_id in known_lnc_set:
                    continue
                if transcript_id in new_lnc_set:
                    lnc_handler.write(str_line)
                elif transcript_id in new_mrna_set:
                    novel_mrnas.add(transcript_id)
                    mrna_handler.write(str_line)
            mrna_ids_handler.write('\n'.join(novel_mrnas))

        return total_lncrans

    def out_fasta(self, identified_new_lnc_list, known_lnc_set, new_mrna_set):
        if not isinstance(identified_new_lnc_list, set):
            identified_new_lnc_list = {i for i in identified_new_lnc_list}

        novel_lnc_file = os.path.join(self.output_dir, 'novel_lncrna.fa')
        novel_mrma_file = os.path.join(self.output_dir, 'novel_mrna.fa')

        with self.option('new_fasta') as in_handler, \
             open(novel_lnc_file, 'w') as lnc_handler, \
             open(novel_mrma_file, 'w') as mrna_handler:
            for seq_id, seq_seq in in_handler:
                trans_id = str(seq_id)
                if trans_id in known_lnc_set:
                    continue
                elif trans_id in identified_new_lnc_list:
                    lnc_handler.write(">" + str(seq_id) + "\n" + str(seq_seq) + "\n")
                elif trans_id in new_mrna_set:
                    mrna_handler.write(">" + str(seq_id) + "\n" + str(seq_seq) + "\n")

    def run(self):
        super(MergePredictionsTool, self).run()
        tools = self.option('tools').strip().split(',')
        preds_dir = self.option('predictions_dir')

        known_lnc_set = {i for i in self.option('known_lnc_json').json_reader()}
        biomart_dict = self.option('biomart_json').json_reader()
        classify_dic = self.classify()

        res_dic = {}  # lncrna detail dict
        extr_list = ['transcript_id', 'score', 'label']
        pfam_extr = ['transcript_id', 'label']
        lncrna_stat = {}  # lncrna statistics in every tool
        trans2gene_info = self.trans2gene_info()
        pfam_predict = None
        for tool in tools:
            file = os.path.join(preds_dir, tool + '_output.txt')
            noncoding_lnc_list, noncoding_gene_list = self.file_reader(
                file, tool, res_dic, trans2gene_info, biomart=biomart_dict, type_dict=classify_dic,
                extr_fields=extr_list if tool != 'pfam' else pfam_extr)
            if tool == 'pfam':
                pfam_predict = file
            noncoding_gene_list = [i for i in set(noncoding_gene_list)]
            noncoding_lnc_list = [i for i in set(noncoding_lnc_list)]
            lncrna_stat[tool] = {
                'tool_name': tool,
                'gene_num': len(noncoding_gene_list),
                'gene_list': noncoding_gene_list,
                'new_lncrna_num': len(noncoding_lnc_list),
                'new_lncrna_list': noncoding_lnc_list
            }
        new_mrna = self.get_new_mrna(pfam_predict)
        # print "new lnc list is {}".format(res_dic)
        identified_new_lnc_list = self.out_put(lncrna_detail_dic=res_dic, lnc_stat_dic=lncrna_stat, tools=tools)
        self.out_new_gtf(identified_new_lnc_list, known_lnc_set=known_lnc_set, new_mrna_set=new_mrna)
        self.out_fasta(identified_new_lnc_list, known_lnc_set=known_lnc_set, new_mrna_set=new_mrna)
        self.end()


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """

        def test(self):
            """
            {'name': 'predictions_dir', 'type': 'string'},
            {'name': 'new_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
            {'name': 'new_fasta', 'type': 'infile', 'format': 'lnc_rna.lnc_fasta'},
            {'name': 'biomart_json', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'classify_info', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'known_lnc_json', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'transcript2gene_info', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'tools', 'type': 'string'},
            {'name': 'identify_num', 'type': 'int', 'default': 1}

            :return:
            """
            import random
            from mbio.workflows.single import SingleWorkflow
            from biocluster.wsheet import Sheet
            data = {
                "id": "MergePredictions_" + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "lnc_rna.lncrna_identification.merge_predictions",
                "instant": False,
                "options": dict(
                    new_gtf="/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/assemble/"
                            "method-cufflinks/output/NewTranscripts/new_transcripts.gtf",
                    new_fasta='/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/assemble/'
                              'method-cufflinks/output/NewTranscripts/new_transcripts.fa',
                    biomart_json='/mnt/ilustre/users/sanger-dev/workspace/20190410/Single_biomart_2370_4297/Biomart/output/biomart.json',
                    classify_info='/mnt/ilustre/users/sanger-dev/workspace/20190410/Single_lncrna_classify_4466_9612/LncrnaClassify/output/novel_lncRNA_classifications.xls',
                    predictions_dir="/mnt/ilustre/users/sanger-dev/workspace/20190312/Single_lncrna_identify_7888_1085/LncrnaIdentify/NewLncrnaPredict/output",
                    # ids_json="/mnt/ilustre/users/sanger-dev/workspace/20190301/Single_BasicFilter_9220/BasicFilter/output/transcript2gene_dict.json",
                    known_lnc_json='/mnt/ilustre/users/sanger-dev/workspace/20190410/Single_extr_known_lnc_9325_7266/ExtrKnownLnc/output/known_lnc_in_new.json',
                    transcript2gene_info="/mnt/ilustre/users/sanger-dev/workspace/20190410/Single_gtf_statistics_8791_868/GtfStat/output/gtf_statistics.json",
                    tools='cpc,cnci,cpat,pfam',
                    identify_num=2,
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()


    unittest.main()
