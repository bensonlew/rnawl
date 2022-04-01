#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/2/27 15:42
@file    : new_lncrna_predict.py
"""
import datetime
import glob

from biocluster.module import Module
import os
import types
from biocluster.core.exceptions import OptionError
import unittest


class NewLncrnaPredictModule(Module):
    def __init__(self, work_id):
        """
        输出文件：
            novel_lncrna.fa
            novel_lncrna.gtf
            novel_lncrna_ids.list
            novel_lncrna_predict_detail.xls
            novel_lncrna_stat.json
            novel_mrna.fa
            novel_mrna.gtf
            novel_mrna_ids.list
        :param work_id:
        """
        super(NewLncrnaPredictModule, self).__init__(work_id)
        options = [
            {'name': 'biomart', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
            {'name': 'biomart_type', 'type': 'string', 'default': 'type1'},
            {'name': 'mrna_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'lnc_db_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},

            {"name": "cpc", "type": "bool", "default": False},
            {"name": "cnci", "type": "bool", "default": False},
            {"name": "cpat", "type": "bool", "default": False},
            {"name": "pfamscan", "type": "bool", "default": False},

            # belong to all tools params
            {'name': 'new_fasta', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'new_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
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
        ]
        self.add_option(options)
        self.step.add_steps('biomart')
        self.biomart_tool = self.add_tool("lnc_rna.lncrna_identification.biomart")
        self.step.add_steps('lnc_classify')
        self.lnc_classify_tool = self.add_tool("lnc_rna.lncrna_identification.lncrna_classify")
        self.step.add_steps('extr_known')
        self.extr_known_tool = self.add_tool("lnc_rna.lncrna_identification.extr_known_lnc")
        self.step.add_steps('basic_filter')
        self.basic_tool = self.add_tool("lnc_rna.lncrna_identification.basic_filter")
        self.step.add_steps('pfam_mrna_filter')
        self.pfam_mrna_filter_tool = self.add_tool("lnc_rna.lncrna_identification.pfam_scan_predictor")
        self.basic_outfile_name = 'basic_filter.fa'
        self.basic_filter_fa = os.path.join(self.basic_tool.output_dir, self.basic_outfile_name)
        self.cpc_tool = None
        self.cpat_tool = None
        self.cnci_tool = None
        self.pfam_tool = None
        self.merge_tool = self.add_tool("lnc_rna.lncrna_identification.merge_predictions")
        self.tools_names = []
        self.predict_tools = []
        self.predict_funcs = []
        self._end_info = 0

    def check_options(self):
        soft_num = 0
        if self.option('cpc') is True:
            soft_num += 1
            self.cpc_tool = self.add_tool("lnc_rna.lncrna_identification.cpc_predictor")
            self.tools_names.append('cpc')
            self.step.add_steps('cpc_predict')
            self.predict_tools.append(self.cpc_tool)
            self.predict_funcs.append(self.cpc_predict)
        if self.option('cnci') is True:
            soft_num += 1
            self.cnci_tool = self.add_tool("lnc_rna.lncrna_identification.cnci_predictor")
            self.tools_names.append('cnci')
            self.step.add_steps('cnci_predict')
            self.predict_tools.append(self.cnci_tool)
            self.predict_funcs.append(self.cnci_predict)
        if self.option('cpat') is True:
            soft_num += 1
            if not (self.option('hexamer_dat').is_set and self.option('logit_model').is_set):
                raise OptionError("使用cpat预测，hexamer_dat 和 logit_model 文件", code="23700801")
            self.cpat_tool = self.add_tool("lnc_rna.lncrna_identification.cpat_predictor")
            self.tools_names.append('cpat')
            self.step.add_steps('cpat_predict')
            self.predict_tools.append(self.cpat_tool)
            self.predict_funcs.append(self.cpat_predict)
        if self.option('pfamscan') is True:
            soft_num += 1
            self.pfam_tool = self.add_tool("lnc_rna.lncrna_identification.pfam_scan_predictor")
            self.step.add_steps('pfamscan_predict')
            self.tools_names.append('pfam')
            self.predict_tools.append(self.pfam_tool)
            self.predict_funcs.append(self.pfamscan_predict)

        # 参数与所选软件不合适则修改
        if soft_num < self.option("identify_num"):
            self.option("identify_num", soft_num)
        self.step.add_steps('predictions_merge')
        if not self.option('new_fasta').is_set:
            raise OptionError("必须输入new_fasta 文件", code="23700801")
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def biomart(self):
        # 输出: biomart.xls, biomart.json
        self.biomart_tool.set_options({
            'biomart': self.option('biomart'),
            'biomart_type': self.option('biomart_type')
        })
        self.biomart_tool.on('start', self.set_step, {'start': self.step.biomart})
        self.biomart_tool.on('end', self.set_step, {'end': self.step.biomart})
        self.biomart_tool.run()

    def classify(self):
        """
            {'type': 'infile', 'name': 'mrna_gtf', 'format': 'lnc_rna.gtf'},
            {'type': 'infile', 'name': 'lncrna_gtf', 'format': 'lnc_rna.gtf'},
            {'type': 'string', 'name': 'out_file_name', 'default': 'lncRNA_classifications.xls'}
        """
        self.lnc_classify_tool.set_options({
            'mrna_gtf': self.option('mrna_gtf').path,
            'lncrna_gtf': self.option('new_gtf').path,
            'out_file_name': 'all_new_rna_classifications.xls'
        })
        self.lnc_classify_tool.on('start', self.set_step, {'start': self.step.lnc_classify})
        self.lnc_classify_tool.on('end', self.set_step, {'end': self.step.lnc_classify})
        self.lnc_classify_tool.run()

    def extr_known(self):
        """
        {'name': 'new_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
        {'name': 'lnc_db_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
        :return: known_lnc_in_new.json
        """
        self.extr_known_tool.set_options({
            'lnc_db_gtf': self.option('lnc_db_gtf').path,
            'new_gtf': self.option('new_gtf').path,
        })
        self.extr_known_tool.on('start', self.set_step, {'start': self.step.extr_known})
        self.extr_known_tool.on('start', self.biomart)
        self.extr_known_tool.on('start', self.classify)
        self.extr_known_tool.on('end', self.set_step, {'end': self.step.extr_known})
        self.extr_known_tool.on('end', self.basic_filter)
        self.extr_known_tool.run()

    def basic_filter(self):
        """
            {'name': 'fasta_file', 'type': 'string'},
            {'name': 'gtf_file', 'type': 'string'},
            {'name': 'transcript_len', 'type': 'int', 'default': 200},
            {'name': 'exon_num', 'type': 'int', 'default': 2},
            {'name': 'orf_len', 'type': 'int', 'default': 300},
        :return:
        """
        self.basic_tool.set_options({
            'fasta_file': self.option('new_fasta').path,
            'gtf_file': self.option('new_gtf').path,
            'transcript_len': self.option('transcript_len'),
            'exon_num': self.option('exon_num'),
            'orf_len': self.option('orf_len'),
            'out_file': self.basic_outfile_name,
            'known_transcripts_json': os.path.join(self.extr_known_tool.output_dir, 'known_lnc_in_new.json'),
        })
        # self.basic_tool.on('end', self.set_output, 'map')
        self.basic_tool.on('start', self.set_step, {'start': self.step.basic_filter})
        self.basic_tool.on('end', self.set_step, {'end': self.step.basic_filter})
        for predict_func in self.predict_funcs:
            self.basic_tool.on('end', predict_func)
        self.basic_tool.on('end', self.mrna_pfam)
        self.basic_tool.run()

    def cpc_predict(self):
        """
            {'name': 'fasta_file', 'type': 'string'},
            {'name': 'cpc_score', 'type': 'float', 'default': 0.5},
        :return:
        """
        self.basic_filter_fa = os.path.join(self.basic_tool.output_dir, self.basic_outfile_name)
        self.cpc_tool.set_options({
            "fasta_file": self.basic_filter_fa,
            "cpc_score": self.option("cpc_score")
        })
        self.cpc_tool.on('start', self.set_step, {'start': self.step.cpc_predict})
        self.cpc_tool.on('end', self.set_output, {'file_regex': 'cpc_output.txt', 'tool_name': 'cpc'})
        self.cpc_tool.on('end', self.set_step, {'end': self.step.cpc_predict})
        self.cpc_tool.run()

    def cnci_predict(self):
        """
            {'name': 'taxonmy', 'type': 'string', 'default': 'Animal'},
            {'name': 'fasta_file', 'type': 'string'},
            {'name': 'cnci_score', 'type': 'float', 'default': 0},
        :return:
        """
        self.basic_filter_fa = os.path.join(self.basic_tool.output_dir, self.basic_outfile_name)
        self.cnci_tool.set_options({
            "taxonmy": self.option('taxonmy'),
            "fasta_file": self.basic_filter_fa,
            "cnci_score": self.option('cnci_score')
        })
        self.cnci_tool.on('start', self.set_step, {'start': self.step.cnci_predict})
        self.cnci_tool.on('end', self.set_output, {'file_regex': 'cnci_output.txt', 'tool_name': 'cnci'})
        self.cnci_tool.on('end', self.set_step, {'end': self.step.cnci_predict})
        self.cnci_tool.run()

    def cpat_predict(self):
        """
            {'name': 'fasta_file', 'type': 'string'},
            {'name': 'hexamer_dat', 'type': 'string', 'default': ''},
            {'name': 'logit_model', 'type': 'string', 'default': ''},
            {'name': 'cpat_score', 'type': 'float', 'default': 0.5}
        :return:
        """
        self.basic_filter_fa = os.path.join(self.basic_tool.output_dir, self.basic_outfile_name)
        self.cpat_tool.set_options({
            "fasta_file": self.basic_filter_fa,
            "hexamer_dat": self.option("hexamer_dat").path,
            "logit_model": self.option("logit_model").path,
            "cpat_score": self.option("cpat_score"),
        })
        self.cpat_tool.on('start', self.set_step, {'start': self.step.cpat_predict})
        self.cpat_tool.on('end', self.set_output, {'file_regex': 'cpat_output.txt', 'tool_name': 'cpat'})
        self.cpat_tool.on('end', self.set_step, {'end': self.step.cpat_predict})
        self.cpat_tool.run()

    def pfamscan_predict(self):
        """
            {'name': 'fasta_file', 'type': 'string'}
        :return:
        """
        self.basic_filter_fa = os.path.join(self.basic_tool.output_dir, self.basic_outfile_name)
        self.pfam_tool.set_options({
            "fasta_file": self.basic_filter_fa,
            "cpu": 12
        })
        self.pfam_tool.on('start', self.set_step, {'start': self.step.pfamscan_predict})
        self.pfam_tool.on('end', self.set_output, {'file_regex': 'pfam_output.txt', 'tool_name': 'pfam'})
        self.pfam_tool.on('end', self.set_step, {'end': self.step.pfamscan_predict})
        self.pfam_tool.run()

    def mrna_pfam(self):
        raw_mrna_fasta_file = os.path.join(self.basic_tool.output_dir, 'raw_mrna.fa')
        self.pfam_mrna_filter_tool.set_options({
            "fasta_file": raw_mrna_fasta_file,
            "cpu": 31
        })
        self.pfam_mrna_filter_tool.on('start', self.set_step, {'start': self.step.pfam_mrna_filter})
        self.pfam_mrna_filter_tool.on('end', self.set_step, {'end': self.step.pfam_mrna_filter})
        self.pfam_mrna_filter_tool.run()

    def predictions_merge(self):
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
        self.merge_tool.set_options({
            'predictions_dir': self.output_dir,
            'tools': ','.join(self.tools_names),
            'new_gtf': self.option('new_gtf').path,
            'new_fasta': self.option('new_fasta').path,
            'identify_num': self.option('identify_num'),
            'biomart_json': os.path.join(self.biomart_tool.output_dir, 'biomart.json'),
            # 'seq_info': os.path.join(self.gtf_stat_tool.output_dir, 'gtf_statistics.json'),
            'classify_info': os.path.join(self.lnc_classify_tool.output_dir, 'all_new_rna_classifications.xls'),
            'known_lnc_json': os.path.join(self.extr_known_tool.output_dir, 'known_lnc_in_new.json'),
            'transcript2gene_info': os.path.join(self.basic_tool.output_dir, 'transcript2gene_dict.json'),
            'pfam_mrna': os.path.join(self.pfam_mrna_filter_tool.output_dir, 'pfam_output.txt')
        })
        self.merge_tool.on('start', self.set_step, {'start': self.step.predictions_merge})
        # self.merge_tool.on('end', self.set_output, {'file_regex': 'new_lncrna*', 'tool_name': 'merge'})
        self.merge_tool.on('end', self.set_output, {'file_regex': 'novel_*', 'tool_name': 'merge'})
        self.merge_tool.on('end', self.set_step, {'end': self.step.predictions_merge})
        self.merge_tool.run()

    def set_output(self, event):
        obj = event['bind_object']
        file_regex = os.path.join(obj.output_dir, event['data']['file_regex'])
        # self.logger.debug('%s ==== file name' % file_regex)
        tool_name = event['data']['tool_name']
        files = glob.glob(file_regex)
        for old_path in files:
            file_name = os.path.basename(old_path)
            new_path = os.path.join(self.output_dir, file_name)
            if os.path.exists(new_path):
                os.remove(new_path)
            os.system('ln {old} {new}'.format(old=old_path, new=new_path))
        self.logger.debug(tool_name + ' has completed prediction')
        if tool_name == 'merge':
            self.end()

    def run(self):
        super(NewLncrnaPredictModule, self).run()
        # for func in self.predict_funcs:
        #     func()
        rely_tools = [self.biomart_tool, self.lnc_classify_tool, self.basic_tool,
                      self.pfam_mrna_filter_tool] + self.predict_tools
        self.on_rely(rely_tools, self.predictions_merge)
        # _ = [self.on_rely(self.basic_tool, tool) for tool in self.predict_tools]
        self.extr_known()


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
                "id": "new_lncrna_predict" + str(random.randint(1, 10000)),
                "type": "module",
                "name": "lnc_rna.new_lncrna_predict",
                "instant": False,
                "options": dict(
                    # {'name': 'biomart', 'type': 'infile', 'format': 'lnc_rna.lnc_common'}
                    biomart='/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/'
                            'Homo_sapiens/Ensemble_release_89/biomart/Homo_sapiens.GRCh38.biomart_gene.txt',
                    #  {'name': 'biomart_type', 'type': 'string', 'default': 'type1'}
                    biomart_type='type1',
                    mrna_gtf='/mnt/ilustre/users/sanger-dev/workspace/20190312/Single_lncrna_identify_7888_1085/'
                             'LncrnaIdentify/KnownLncrnaIdentify/output/mrna.gtf',
                    # {'name': 'lnc_db_gtf', 'type': 'infile', 'format': 'lnc_rna.lnc_common'},
                    lnc_db_gtf='/mnt/ilustre/users/sanger-dev/workspace/20190409/LncDb_lnc_db_workflow_3956_9130/'
                               'output/lncrna.gtf',
                    # {"name": "cpc", "type": "bool", "default": False}
                    cpc=True,
                    # {"name": "cnci", "type": "bool", "default": False}
                    cnci=True,
                    # {"name": "cpat", "type": "bool", "default": False}
                    cpat=True,
                    # {"name": "pfamscan", "type": "bool", "default": False}
                    pfamscan=True,
                    # {'name': 'new_fasta', 'type': 'infile', 'format': 'lnc_rna.fasta'},
                    new_fasta=r'/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/assemble/'
                              r'method-cufflinks/output/NewTranscripts/new_transcripts.fa',
                    # {'name': 'new_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
                    new_gtf="/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/assemble/"
                            "method-cufflinks/output/NewTranscripts/new_transcripts.gtf",
                    # {'name': 'identify_num', 'type': 'int', 'default': 2}
                    identify_num=2,
                    # {'name': 'transcript_len', 'type': 'int', 'default': 200}
                    transcript_len=200,
                    # {'name': 'exon_num', 'type': 'int', 'default': 2}
                    exon_num=2,
                    # {'name': 'orf_len', 'type': 'int', 'default': 300}
                    orf_len=300,
                    # {'name': 'cpc_score', 'type': 'float', 'default': 0.5}
                    cpc_score=0.5,
                    # {'name': 'taxonmy', 'type': 'string', 'default': 'Animal'}
                    taxonmy='Human',
                    # {'name': 'cnci_score', 'type': 'float', 'default': 0}
                    cnci_score=0,
                    # {'name': 'cpat_score', 'type': 'float', 'default': 0.5}
                    cpat_score=0.5,
                    # {'name': 'hexamer_dat', 'type': 'infile', 'format': 'lnc_rna.common'}
                    hexamer_dat="/mnt/ilustre/users/isanger/app/bioinfo/lnc_rna/CPAT-1.2.4/dat/Human_Hexamer.tsv",
                    # {'name': 'logit_model', 'type': 'infile', 'format': 'lnc_rna.common'}
                    logit_model="/mnt/ilustre/users/isanger/app/bioinfo/lnc_rna/CPAT-1.2.4/dat/Human_logitModel.RData",
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()


    unittest.main()
