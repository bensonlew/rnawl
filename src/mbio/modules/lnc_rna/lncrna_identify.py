#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/3/7 22:52
@file    : lncrna_identify.py
"""

import os
import unittest
from itertools import chain
from pathlib import Path

from biocluster.core.exceptions import OptionError, FileError
from biocluster.module import Module

from mbio.files.lnc_rna.lnc_fasta import LncFastaFile
from mbio.files.lnc_rna.lnc_gtf import LncGtfFile


class LncrnaIdentifyModule(Module):
    def __init__(self, work_id):
        super(LncrnaIdentifyModule, self).__init__(work_id)
        options = [
            # new lncrna prediction params
            {"name": "cpc", "type": "bool", "default": False},
            {"name": "cnci", "type": "bool", "default": False},
            {"name": "cpat", "type": "bool", "default": False},
            {"name": "pfamscan", "type": "bool", "default": False},
            # belong to all tools params
            {'name': 'assemble_dir', 'type': 'string', 'required': True},
            {'name': 'identify_num', 'type': 'int', 'default': 2},
            # basic filter tool params
            {'name': 'transcript_len', 'type': 'int', 'default': 200},
            {'name': 'exon_num', 'type': 'int', 'default': 2},
            {'name': 'orf_len', 'type': 'int', 'default': 300},
            # cpc params
            {'name': 'cpc_score', 'type': 'float', 'default': 0.5},
            # cnci params
            {'name': 'cnci_score', 'type': 'float', 'default': 0},
            {'name': 'taxonmy', 'type': 'string', 'default': 'Animal'},
            # cpat parmas
            {'name': 'cpat_score', 'type': 'float', 'default': 0.5},
            {'name': 'hexamer_dat', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'logit_model', 'type': 'infile', 'format': 'lnc_rna.common'},

            {'name': 'rna_exp_matrix', 'type': 'infile', 'format': 'lnc_rna.common'},

            # known lncrna identify
            {'name': 'lnc_db_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
            {'name': 'lnc_ids_file', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'rna_exp_matrix', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'biomart', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'biomart_type', 'type': 'string', 'default': 'type2'},
            {'name': 'lnc_ref_db', 'type': 'string', 'default': 'ensembl'},

            # 是否只保留表达的rna, 默认True
            {'name': 'retain_expressed_mrna', 'type': 'bool', 'default': True},
        ]

        self.add_option(options)
        self.step.add_steps('known_lncrna_identify',
                            'known_lncrna_classify',
                            'novel_lncrna_classify',
                            'lncrna_stat',
                            'new_lncrna_predict')

        self.known_lncrna_tool = self.add_tool("lnc_rna.lncrna_identification.known_lncrna_identify")
        self.known_lnc_out_path = Path(self.known_lncrna_tool.output_dir)

        self.known_lncrna_classify_tool = self.add_tool("lnc_rna.lncrna_identification.lncrna_classify")
        self.known_cla_path = Path(self.known_lncrna_classify_tool.output_dir)

        self.novel_lncrna_classify_tool = self.add_tool("lnc_rna.lncrna_identification.lncrna_classify")
        self.novel_cla_path = Path(self.novel_lncrna_classify_tool.output_dir)

        self.lncrna_stat_tool = self.add_tool("lnc_rna.lncrna_identification.lncrna_stat")
        self.lnc_stat_path = Path(self.lncrna_stat_tool.output_dir)

        self.new_predict_tool = self.add_module("lnc_rna.new_lncrna_predict")
        self.new_predict_path = Path(self.new_predict_tool.output_dir)

        self._end_info = 0
        # 组装后文件夹路径
        self.assenble_outdir_path = None

    def check_options(self):
        data_dir = self.option('assemble_dir')
        check_files = ['all_transcripts.fa', 'new_transcripts.fa', 'ref_and_new.gtf', 'new_transcripts.gtf',
                       'new_genes.gtf']
        for file in (os.path.join(data_dir, n) for n in check_files):
            if not os.path.isfile(file):
                raise FileError('File not found: %s' % file)
        if self.option('cpat') is True:
            if not self.option('hexamer_dat').is_set:
                raise OptionError('hexamer_dat is not set when cpat is True')
            if not self.option('logit_model').is_set:
                raise OptionError('logit_model is not set when cpat is True')

        self.assenble_outdir_path = Path(self.option('assemble_dir'))
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def known_lncrna_filter(self):
        """
            {'name': 'all_fasta', 'type': 'infile', 'format': 'lnc_rna.lnc_fasta'},
            {'name': 'new_fasta', 'type': 'infile', 'format': 'lnc_rna.lnc_fasta'},
            {'name': 'new_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
            {'name': 'ref_and_new_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
            {'name': 'lnc_db_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
            {'name': 'lnc_ids_file', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'rna_exp_matrix', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'biomart', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'lnc_ref_db', 'type': 'string', 'default': 'ensembl'},

            输出：
                mrna.gtf
                known_lncrna.fa
                known_lncrna.gtf
                known_lncrna_ids.list
                known_lncrna_info.txt：统计信息
                new_transcripts_list.json: 过滤掉已知转录本后的new_transcripts_list
                all_transcript2gene_info.json: 转录本对应基因详情 -- gene_id, gene_name, ...
        :return:
        """
        file_objs = LncFastaFile(), LncFastaFile(), LncGtfFile(), LncGtfFile()
        file_names = 'all_transcripts.fa', 'new_transcripts.fa', 'new_transcripts.gtf', 'ref_and_new.gtf'
        for file_obj, name in zip(file_objs, file_names):
            file_obj.set_path(str(self.assenble_outdir_path.joinpath(name)))
        # files = [str(self.assenble_outdir_path.joinpath(name)) for name in file_names]
        # all_fasta, new_fasta, new_gtf, ref_and_new_gtf = files
        all_fasta, new_fasta, new_gtf, ref_and_new_gtf = file_objs
        self.known_lncrna_tool.set_options(dict(
            all_fasta=all_fasta,
            new_fasta=new_fasta,
            new_gtf=new_gtf,
            ref_and_new_gtf=ref_and_new_gtf,
            lnc_db_gtf=self.option('lnc_db_gtf'),
            lnc_ids_file=self.option('lnc_ids_file'),
            rna_exp_matrix=self.option('rna_exp_matrix'),
            biomart=self.option('biomart'),
            biomart_type=self.option('biomart_type'),
            lnc_ref_db=self.option('lnc_ref_db')
        ))
        # self.basic_tool.on('end', self.set_output, 'map')
        self.known_lncrna_tool.on('start', self.set_step, {'start': self.step.known_lncrna_identify})
        self.known_lncrna_tool.on('end', self.set_step, {'end': self.step.known_lncrna_identify})

        # new lncrna 中去除已知后剩余的转录本
        new_lnc_json = self.known_lnc_out_path.joinpath('new_transcripts_list.json')
        # all_transcript2gene_info.json
        transcript2gene_info = self.known_lnc_out_path.joinpath('all_transcript2gene_info.json')
        self.known_lncrna_tool.on('end', self.new_lncrna_predict,
                                  {'new_lnc_json': str(new_lnc_json), 'transcript2gene_info': str(transcript2gene_info)})

        # 已知lncRNA GTF文件路径
        lnc_gtf_path = self.known_lnc_out_path.joinpath('known_lncrna.gtf')
        # 已知lncRNA分类
        self.known_lncrna_tool.on('end', self.lncrna_classify, {'lnc_gtf_path': str(lnc_gtf_path), 'is_novel': False})

        self.known_lncrna_tool.run()

    def new_lncrna_predict(self, event):
        """
            {"name": "cpc", "type": "bool", "default": False},
            {"name": "cnci", "type": "bool", "default": False},
            {"name": "cpat", "type": "bool", "default": False},
            {"name": "pfamscan", "type": "bool", "default": False},

            # belong to all tools params
            {'name': 'fasta_file', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'gtf_file', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'identify_num', 'type': 'int', 'default': 2},

            # basic filter tool params
            {'name': 'transcript_len', 'type': 'int', 'default': 200},
            {'name': 'exon_num', 'type': 'int', 'default': 2},
            {'name': 'orf_len', 'type': 'int', 'default': 300},
            {'name': 'new_transcripts_json_file', 'type': 'string'},

            # cpc params
            {'name': 'cpc_score', 'type': 'float', 'default': 0.5},
            # cnci params
            {'name': 'cnci_score', 'type': 'float', 'default': 0},
            {'name': 'taxonmy', 'type': 'string', 'default': 'Animal'},
            # cpat parmas
            {'name': 'cpat_score', 'type': 'float', 'default': 0.5},
            {'name': 'hexamer_dat', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'logit_model', 'type': 'infile', 'format': 'lnc_rna.common'},
            # {'name': 'out_file', 'type': 'string', 'default': 'basic_filter_1.fa'},
            # merge params
            {'name': 'transcript2gene_info', 'type': 'infile', 'format': 'lnc_rna.common'}

            输出文件：
                new_lncrna.gtf     :  预测到的lncrna GTF文件
                new_lncrna_ids.list           :  新lncrna ids列表
                new_lncrna_stat.json          :  新lncrna统计结果文件 [用于导表]
                new_lncrna_predict_detail.json:  新lncrna预测详情表 [用于导表]
        :return:
        """
        data = event['data']
        self.new_predict_tool.set_options(dict(
            cpc=self.option('cpc'),
            cnci=self.option('cnci'),
            cpat=self.option('cpat'),
            pfamscan=self.option('pfamscan'),
            fasta_file=str(self.assenble_outdir_path.joinpath('new_transcripts.fa')),
            gtf_file=str(self.assenble_outdir_path.joinpath('new_transcripts.gtf')),
            identify_num=self.option('identify_num'),
            transcript_len=self.option('transcript_len'),
            exon_num=self.option('exon_num'),
            orf_len=self.option('orf_len'),
            taxonmy=self.option('taxonmy'),
            cpc_score=self.option('cpc_score'),
            cnci_score=self.option('cnci_score'),
            cpat_score=self.option('cpat_score'),
            hexamer_dat=self.option('hexamer_dat'),
            logit_model=self.option('logit_model'),
            new_transcripts_json_file=data['new_lnc_json'],
            transcript2gene_info=data['transcript2gene_info']
        ))
        self.new_predict_tool.on('start', self.set_step, {'start': self.step.new_lncrna_predict})
        self.new_predict_tool.on('end', self.set_step, {'end': self.step.new_lncrna_predict})

        # 新lncRNA分类, 输出的新lncrna GTF文件路径
        lnc_gtf_path = self.new_predict_path.joinpath('new_lncrna.gtf')
        self.new_predict_tool.on('end', self.lncrna_classify, {'lnc_gtf_path': str(lnc_gtf_path), 'is_novel': True})

        self.new_predict_tool.run()

    def lncrna_classify(self, event):
        """
            {'type': 'infile', 'name': 'mrna_gtf', 'format': 'lnc_rna.gtf'},
            {'type': 'infile', 'name': 'lncrna_gtf', 'format': 'lnc_rna.gtf'},
            {'type': 'string', 'name': 'out_file_name', 'default': 'lncRNA_classifications.xls'}

            输出文件：
                novel_lncRNA_classifications.xls: lncrna分类文件
                known_lncRNA_classifications.xls: lncrna分类文件

        :param lnc_gtf:
        :return:
        """
        data = event['data']
        lnc_gtf_path, is_novel = data['lnc_gtf_path'], data['is_novel']
        mrna_gtf_path = os.path.join(self.known_lncrna_tool.output_dir, 'mrna.gtf')
        self.logger.debug('mrna_gtf_path ' + mrna_gtf_path)
        # mrna_gtf, lnc_gtf = GtfFile(), GtfFile()
        # mrna_gtf.set_path(mrna_gtf_path)
        # lnc_gtf.set_path(lnc_gtf_path)

        self.logger.debug('lnc_gtf_path' + lnc_gtf_path + '---' + str(is_novel))
        if is_novel:
            # novel_lncRNA_classifications.xls: lncrna分类文件
            self.novel_lncrna_classify_tool.set_options(dict(
                mrna_gtf=mrna_gtf_path,
                lncrna_gtf=lnc_gtf_path,
                out_file_name='novel_lncRNA_classifications.xls'
            ))
            self.novel_lncrna_classify_tool.on('start', self.set_step, {'start': self.step.novel_lncrna_classify})
            self.novel_lncrna_classify_tool.on('end', self.set_step, {'end': self.step.novel_lncrna_classify})

            # 运行 lncrna 统计
            # 新lncrna分类结果文件
            novel_classify = self.novel_cla_path.joinpath('novel_lncRNA_classifications.xls')
            # 已知lncrna分类结果文件
            known_classify = self.known_cla_path.joinpath('known_lncRNA_classifications.xls')
            # 新lncrna ids列表
            new_lncrna_ids = self.new_predict_path.joinpath('new_lncrna_ids.list')
            # 已知lncrna ids列表
            known_lncrna_ids = self.known_lnc_out_path.joinpath('known_lncrna_ids.list')
            # 转录本对应基因信息映射
            transcript2gene_info = self.known_lnc_out_path.joinpath('all_transcript2gene_info.json')

            data = {
                'novel_classify': str(novel_classify),
                'known_classify': str(known_classify),
                'new_lncrna_ids': str(new_lncrna_ids),
                'known_lncrna_ids': str(known_lncrna_ids),
                'transcript2gene_info': str(transcript2gene_info)
            }

            self.novel_lncrna_classify_tool.on('end', self.lncrna_stat, data)

            self.novel_lncrna_classify_tool.run()
        else:
            # known_lncRNA_classifications.xls: lncrna分类文件
            self.known_lncrna_classify_tool.set_options(dict(
                mrna_gtf=mrna_gtf_path,
                lncrna_gtf=lnc_gtf_path,
                out_file_name='known_lncRNA_classifications.xls'
            ))
            self.known_lncrna_classify_tool.on('start', self.set_step, {'start': self.step.known_lncrna_classify})
            self.known_lncrna_classify_tool.on('end', self.set_step, {'end': self.step.known_lncrna_classify})
            self.known_lncrna_classify_tool.run()

    def lncrna_stat(self, event):
        """
            {'name': 'new_lncrna_ids', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'known_lncrna_ids', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'rna_exp_matrix', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'novel_classify_file', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'known_classify_file', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'new_genes_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
            {'name': 'ref_new_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
            {'name': 'all_trans_fa', 'type': 'infile', 'format': 'lnc_rna.lnc_fasta'},
            {'name': 'transcript2gene_info', 'type': 'infile', 'format': 'lnc_rna.common'},

            # data 格式
                {'novel_classify': novel_classify, 'known_classify': known_classify,
                'new_lncrna_ids': new_lncrna_ids, 'known_lncrna_ids': known_lncrna_ids}
            # 输出文件：
                novel_lncrna.gtf
                known_lncrna.gtf
                mrna_lncrna.gtf

                novel_lncrna.fa
                known_lncrna.fa
                mrna_lncrna.fa
        :return:
        """
        data = event['data']

        self.lncrna_stat_tool.set_options(dict(
            new_lncrna_ids=data['new_lncrna_ids'],
            known_lncrna_ids=data['known_lncrna_ids'],
            rna_exp_matrix=self.option('rna_exp_matrix'),
            novel_classify_file=data['novel_classify'],
            known_classify_file=data['known_classify'],
            # 对rna是否表达进行过滤
            # {'name': 'total_known_lnc_ids', 'type': 'infile', 'format': 'lnc_rna.common'},
            # {'name': 'expressed_trans_ids', 'type': 'infile', 'format': 'lnc_rna.common'},
            # {'name': 'retain_expressed_mrna', 'type': 'bool', 'default': True},
            total_known_lnc_ids=str(self.known_lnc_out_path.joinpath('total_known_lnc_ids.json')),
            expressed_trans_ids=str(self.known_lnc_out_path.joinpath('expressed_transcripts_ids.json')),
            retain_expressed_mrna=self.option('retain_expressed_mrna'),
            # new lncrna 中去除已知后剩余的转录本
            filtered_new_transcripts=str(self.known_lnc_out_path.joinpath('new_transcripts_list.json')),
            total_new_transcripts=str(self.known_lnc_out_path.joinpath('total_new_transcripts_list.json')),
            new_genes_gtf=str(self.assenble_outdir_path.joinpath('new_genes.gtf')),
            ref_new_gtf=str(self.assenble_outdir_path.joinpath('ref_and_new.gtf')),
            all_trans_fa=str(self.assenble_outdir_path.joinpath('all_transcripts.fa')),
            transcript2gene_info=data['transcript2gene_info']
        ))
        self.lncrna_stat_tool.on('start', self.set_step, {'start': self.step.lncrna_stat})
        self.lncrna_stat_tool.on('end', self.set_step, {'end': self.step.lncrna_stat})
        self.lncrna_stat_tool.on('end', self.set_output)
        self.lncrna_stat_tool.run()

    def set_output(self):
        """
            novel_out_gtf = os.path.join(self.output_dir, 'novel_lncrna.gtf')
            known_out_gtf = os.path.join(self.output_dir, 'known_lncrna.gtf')
            lnc_out_gtf = os.path.join(self.output_dir, 'total_lncrna.gtf')
            novel_mrna_out_gtf = os.path.join(self.output_dir, 'novel_mrna.gtf')
            known_mrna_out_gtf = os.path.join(self.output_dir, 'known_mrna.gtf')
            total_mrna_out_gtf = os.path.join(self.output_dir, 'total_mrna.gtf')

            novel_lnc_out_fa = os.path.join(self.output_dir, 'novel_lncrna.fa')
            known_lnc_out_fa = os.path.join(self.output_dir, 'known_lncrna.fa')
            total_lnc_out_fa = os.path.join(self.output_dir, 'total_lncrna.fa')
            novel_mrna_out_fa = os.path.join(self.output_dir, 'novel_mrna.fa')
            known_mrna_out_fa = os.path.join(self.output_dir, 'known_mrna.fa')
            total_mrna_out_fa = os.path.join(self.output_dir, 'total_mrna.fa')
        :return:
        """
        known_out = ['known_lncrna.fa', 'known_lncrna.gtf', 'known_lncrna_ids.list', 'known_lncrna_info.txt']
        novel_out = ['new_lncrna.gtf', 'new_lncrna_ids.list', 'new_lncrna_stat.json', 'new_lncrna_predict_detail.json']
        lncrna_stat_out = ['lncrna_stat_in_samples.json', 'lncrna_classify_stat.json',
                           # 输出 gtf 文件
                           'novel_lncrna.gtf', 'known_lncrna.gtf', 'total_lncrna.gtf',
                           'novel_mrna.gtf', 'known_mrna.gtf', 'total_mrna.gtf',
                           # 输出 fasta 文件
                           'novel_lncrna.fa', 'known_lncrna.fa', 'total_lncrna.fa',
                           'novel_mrna.fa', 'known_mrna.fa', 'total_mrna.fa',
                           # lnc 和 mrna id列表
                           'total_lncrna.list', 'total_mrna.list']

        out_data = chain(
            ((self.known_lnc_out_path, out) for out in known_out),
            ((self.new_predict_path, out) for out in novel_out),
            ((self.lnc_stat_path, out) for out in lncrna_stat_out),
            ((self.known_cla_path, 'known_lncRNA_classifications.xls'),),
            ((self.novel_cla_path, 'novel_lncRNA_classifications.xls'),)
        )
        for out_path, file_name in out_data:
            f_path = out_path.joinpath(file_name)
            new_path = Path(self.output_dir, f_path.name)

            if new_path.is_file():
                new_path.unlink()

            os.system('ln {old} {new}'.format(old=str(f_path), new=str(new_path)))

        self.end()

    def run(self):
        super(LncrnaIdentifyModule, self).run()
        self.known_lncrna_filter()


if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """

        def test(self):
            """
            {"name": "cpc", "type": "bool", "default": False},
            {"name": "cnci", "type": "bool", "default": False},
            {"name": "cpat", "type": "bool", "default": False},
            {"name": "pfamscan", "type": "bool", "default": False},

            # belong to all tools params
            {'name': 'fasta_file', 'type': 'string', 'format': 'lnc_rna.fasta'},
            {'name': 'gtf_file', 'type': 'string', 'format': 'lnc_rna.gtf'},

            # basic filter tool params
            {'name': 'transcript_len', 'type': 'int', 'default': 200},
            {'name': 'exon_num', 'type': 'int', 'default': 2},
            {'name': 'orf_len', 'type': 'int', 'default': 300},

            # cpc params
            {'name': 'cpc_score', 'type': 'float', 'default': 0.5},
            # cnci params
            {'name': 'cnci_score', 'type': 'float', 'default': 0},
            {'name': 'taxonmy', 'type': 'string', 'default': 'Animal'},
            # cpat parmas
            {'name': 'cpat_score', 'type': 'float', 'default': 0.5},
            {'name': 'hexamer_dat', 'type': 'string', 'default': 'lnc_rna.common'},
            {'name': 'logit_model', 'type': 'string', 'default': 'lnc_rna.common'},
            :return:
            """
            import random
            from mbio.workflows.single import SingleWorkflow
            from biocluster.wsheet import Sheet
            data = {
                "id": "lncrna_identify_" + str(random.randint(1, 10000)) + '_' + str(random.randint(1, 10000)),
                "type": "module",
                "name": "lnc_rna.lncrna_identify",
                "instant": False,
                "options": dict(
                    # {"name": "cpc", "type": "bool", "default": False},
                    # {"name": "cnci", "type": "bool", "default": False},
                    # {"name": "cpat", "type": "bool", "default": False},
                    # {"name": "pfamscan", "type": "bool", "default": False},
                    cpc=True,
                    cnci=True,
                    cpat=True,
                    pfamscan=True,
                    # belong to all tools params
                    # {'name': 'assemble_dir', 'type': 'string'},
                    assemble_dir='/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/assemble/method-cufflinks/output/NewTranscripts/',
                    # {'name': 'identify_num', 'type': 'int', 'default': 2},
                    identify_num=2,
                    # basic filter tool params
                    # {'name': 'transcript_len', 'type': 'int', 'default': 200},
                    transcript_len=200,
                    # {'name': 'exon_num', 'type': 'int', 'default': 2},
                    exon_num=2,
                    # {'name': 'orf_len', 'type': 'int', 'default': 300},
                    orf_len=300,
                    # cpc params
                    # {'name': 'cpc_score', 'type': 'float', 'default': 0.5},
                    cpc_score=0.5,
                    # cnci params
                    # {'name': 'cnci_score', 'type': 'float', 'default': 0},
                    cnci_score=0,
                    # {'name': 'taxonmy', 'type': 'string', 'default': 'Animal'},
                    taxonmy='Human',
                    # cpat parmas
                    # {'name': 'cpat_score', 'type': 'float', 'default': 0.5},
                    cpat_score=0.5,
                    # {'name': 'hexamer_dat', 'type': 'infile', 'format': 'lnc_rna.common'},
                    # {'name': 'logit_model', 'type': 'infile', 'format': 'lnc_rna.common'},
                    hexamer_dat="/mnt/ilustre/users/isanger/app/bioinfo/lnc_rna/CPAT-1.2.4/dat/Human_Hexamer.tsv",
                    logit_model="/mnt/ilustre/users/isanger/app/bioinfo/lnc_rna/CPAT-1.2.4/dat/Human_logitModel.RData",
                    # {'name': 'rna_exp_matrix', 'type': 'infile', 'format': 'lnc_rna.common'},
                    rna_exp_matrix="/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncrna_tools/transcript.tpm.matrix",

                    # known lncrna identify
                    # {'name': 'lnc_db_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_gtf'},
                    # {'name': 'lnc_ids_file', 'type': 'infile', 'format': 'lnc_rna.common'},
                    # {'name': 'biomart', 'type': 'infile', 'format': 'lnc_rna.common'},
                    # {'name': 'lnc_ref_db', 'type': 'string', 'default': 'ensembl'}
                    lnc_db_gtf="/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/human/38/ensembl/lnc_db/Homo_sapiens.GRCh38.95.gtf",
                    lnc_ids_file="/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/lncRNA_2019_01_09/human/38/ensembl/lnc_db/ids_matrix.xls",
                    biomart=r'/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/biomart/Homo_sapiens.GRCh38.biomart_gene.txt',
                    biomart_type='type1',
                    lnc_ref_db='ensembl',

                    retain_expressed_mrna=True,
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()


    unittest.main()
