#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/2/26 15:35
@file    : cpat_predictor.py
"""
import csv
import os
import unittest

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
from mbio.packages.lnc_rna.lnc_identification.predict_common_tools import CommandFactory


class CpatPredictorAgent(Agent):
    def __init__(self, parent):
        super(CpatPredictorAgent, self).__init__(parent)
        options = [
            {'name': 'fasta_file', 'type': 'string'},
            {'name': 'hexamer_dat', 'type': 'string', 'default': ''},
            {'name': 'logit_model', 'type': 'string', 'default': ''},
            {'name': 'cpat_score', 'type': 'float', 'default': 0.5},
            {'name': 'cpu', 'type': 'int', 'default': 2}
        ]
        cpu = 10
        self.add_option(options)
        self.option('cpu', cpu)
        self.step.add_steps("cpat_predict")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.cpat_predict.start()
        self.step.update()

    def step_end(self):
        self.step.cpat_predict.finish()
        self.step.update()

    def check_options(self):
        for name in ('fasta_file', 'hexamer_dat', 'logit_model'):
            file = self.option(name)
            if file and not os.path.isfile(file):
                raise OptionError('缺少fasta file文件')
        return True

    def set_resource(self):
        self._cpu = 3
        self._memory = '5G'

    def end(self):
        # self.add_upload_dir(self.output_dir)
        super(CpatPredictorAgent, self).end()


class CpatPredictorTool(Tool):
    def __init__(self, config):
        super(CpatPredictorTool, self).__init__(config)
        self.python_path = "miniconda2/bin/python"
        env_path = ':'.join([
            os.path.join(self.config.SOFTWARE_DIR, 'program/R-3.3.3/bin/'),
            os.path.join(self.config.SOFTWARE_DIR, 'program/Python/bin')
        ])
        self.set_environ(PATH=env_path)
        self.cpat_out_file = os.path.join(self.output_dir, 'cpat_output.txt')

    def cmd_runner(self, cmd_name, cmd, check_stat=True):
        cmd_obj = self.add_command(cmd_name, cmd)
        cmd_obj.run()
        self.wait()
        if check_stat is True:
            self.check_stat(cmd_obj)
        return cmd_obj

    def check_stat(self, *cmd_objs):
        self.wait(*cmd_objs)
        for cmd_obj in cmd_objs:
            if cmd_obj.return_code == 0:
                self.logger.info('%s：运行完成' % cmd_obj.cmd)
            elif cmd_obj.return_code in (1, -9):
                self.add_state('memory_limit', 'memory is low!')
            else:
                self.set_error('%s: 运行错误%s' % cmd_obj.cmd)

    def cpat_predict(self):
        """
        -g GENE_FILE, --gene=GENE_FILE
                        Transcripts either in BED format or mRNA sequences in
                        FASTA format: If this is BED format file, '-r' must be
                        specified; if this is mRNA sequence file in FASTA
                        format, ignore the '-r' option. The input BED or FASTA
                        file could be regular text file or compressed file
                        (*.gz, *.bz2) or accessible url.
        -o OUT_FILE, --outfile=OUT_FILE
                            output file. Tab separated text file: geneID <tab>
                            mRNA size <tab> ORF size <tab> Fickett Score <tab>
                            Hexamer Score<tab>Coding Probability.
        -x HEXAMER_DAT, --hex=HEXAMER_DAT
                            Prebuilt hexamer frequency table (Human, Mouse, Fly,
                            Zebrafish). Run 'make_hexamer_tab.py' to make this
                            table out of your own training dataset.
        -d LOGIT_MODEL, --logitModel=LOGIT_MODEL
                            Prebuilt training model (Human, Mouse, Fly,
                            Zebrafish). Run 'make_logitModel.py' to build logit
                            model out of your own training datset
        -r REF_GENOME, --ref=REF_GENOME
                            Reference genome sequences in FASTA format. Ignore
                            this option if mRNA sequences file was provided to
                            '-g'. Reference genome file will be indexed
                            automatically (produce *.fai file along with the
                            original *.fa file within the same directory) if
                            hasn't been done.
        -s START_CODONS, --start=START_CODONS
                            Start codon (DNA sequence, so use 'T' instead of 'U')
                            used to define open reading frame (ORF). default=ATG
        -t STOP_CODONS, --stop=STOP_CODONS
                            Stop codon (DNA sequence, so use 'T' instead of 'U')
                            used to define open reading frame (ORF). Multiple stop
                            codons should be separated by ','. default=TAG,TAA,TGA
        :return:
        """
        # '{soft} -g {fa_file} -d {logit_model} -x {hexamer_dat} -o {outfile}'
        cnci_tool = os.path.join(self.config.SOFTWARE_DIR, 'program/Python/bin/cpat.py')
        cpat_out_file = os.path.join(self.work_dir, 'temp_cpat_output.txt')

        cmd_fac_obj = CommandFactory(self.python_path + ' ' + cnci_tool)
        cmd_fac_obj.add_param('-g', self.option('fasta_file'), param_desc='fasta file path')
        # logit model example: Human_logitModel.RData
        cmd_fac_obj.add_param('-d', self.option('logit_model'), param_desc='logit model')
        # hexamer dat, example: Human_Hexamer.tsv
        cmd_fac_obj.add_param('-x', self.option('hexamer_dat'), param_desc='hexamer dat')
        cmd_fac_obj.add_param('-o', cpat_out_file, param_desc='output file path')

        self.cmd_runner('cpat_predict', cmd_fac_obj.to_string(), check_stat=True)

        return cpat_out_file

    def lncrna_mark(self, predict_file):
        # mRNA_size  ORF_size  Fickett_score  Hexamer_score  coding_prob
        out_fields = ['transcript_id', 'transcript_length', 'orf_length', 'score', 'label']
        extr_fields = ['transcript_id', 'mRNA_size', 'ORF_size', 'score', 'label']
        out_demo = '\t'.join('{%s}' % i for i in extr_fields) + '\n'
        with open(predict_file) as in_handler, open(self.cpat_out_file, 'w') as out_handler:
            fields = in_handler.readline().strip().split('\t')
            fields.insert(0, 'transcript_id')
            out_handler.write('\t'.join(out_fields) + '\n')
            for line_list in csv.reader(in_handler, delimiter='\t'):
                line_dic = {k: v for k, v in zip(fields, line_list)}
                coding_prob = float(line_dic['coding_prob'])
                line_dic['label'] = 'noncoding' if coding_prob < self.option('cpat_score') else 'coding'
                out_handler.write(out_demo.format(score=coding_prob, **line_dic))

    def run(self):
        super(CpatPredictorTool, self).run()
        cpat_out_file = self.cpat_predict()
        self.lncrna_mark(cpat_out_file)
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
                "id": "CpatPredictor_" + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "tool_lab.coding_prediction.cpat_predictor",
                "instant": False,
                "options": dict(
                    fasta_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab_batch2/Coding_predict/basic_filter.fa",
                    hexamer_dat="/mnt/ilustre/users/isanger/app/bioinfo/lnc_rna/CPAT-1.2.4/dat/Human_Hexamer.tsv",
                    logit_model="/mnt/ilustre/users/isanger/app/bioinfo/lnc_rna/CPAT-1.2.4/dat/Human_logitModel.RData",
                    cpat_score=0.5,
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()

    unittest.main()
