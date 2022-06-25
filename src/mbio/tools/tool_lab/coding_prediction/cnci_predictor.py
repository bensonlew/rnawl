#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/2/26 14:28
@file    : cnci_predictor.py
"""
import csv
import os
import unittest

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
from mbio.packages.lnc_rna.lnc_identification.predict_common_tools import CommandFactory


class CnciPredictorAgent(Agent):
    def __init__(self, parent):
        super(CnciPredictorAgent, self).__init__(parent)
        # 'fasta_file', 'gtf_file', 'hexamer_dat', 'logit_model',
        # , 'transcript_len', 'cnci_score',
        # 'orf_len', 'cpat_score', 'taxonmy', 'exon_num'
        options = [
            {'name': 'taxonmy', 'type': 'string', 'default': 'Animal'},
            {'name': 'fasta_file', 'type': 'string'},
            {'name': 'cnci_score', 'type': 'float', 'default': 0},
            {'name': 'cpu', 'type': 'int', 'default': 2}
        ]
        self.add_option(options)
        self.step.add_steps("cnci_predict")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.cnci_predict.start()
        self.step.update()

    def step_end(self):
        self.step.cnci_predict.finish()
        self.step.update()

    def check_options(self):
        fa_file = self.option('fasta_file')
        if fa_file and not os.path.isfile(fa_file):
            raise OptionError('缺少fasta file文件')
        return True

    def set_resource(self):
        self._cpu = 12
        self._memory = '8G'
        self.option('cpu', self._cpu)

    def end(self):
        # self.add_upload_dir(self.output_dir)
        super(CnciPredictorAgent, self).end()


class CnciPredictorTool(Tool):
    def __init__(self, config):
        super(CnciPredictorTool, self).__init__(config)
        self.python_path = "miniconda2/bin/python"
        self.env_path = self.config.SOFTWARE_DIR + '/miniconda2/bin'
        self.set_environ(PATH=self.env_path)
        self.cnci_out_file = os.path.join(self.output_dir, 'cnci_output.txt')

    def cmd_runner(self, cmd_name, cmd, check_stat=True):
        cmd_obj = self.add_command(cmd_name, cmd)
        cmd_obj.run()
        # self.wait(cmd_obj)
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

    def cnci_predict(self):
        """
         -f input files, --file=input files
                        enter your transcript (sequence or gtf)
          -o output files, --out=output files
                                assign your output file
          -p prallel numbers, --parallel=prallel numbers
                                please enter your specified speed ratio
          -m model types, --model=model types
        :return:
        """
        taxonmy = self.option('taxonmy')
        model = 'pl' if taxonmy.lower() == 'plant' else 've'
        # software_dir = /mnt/ilustre/users/sanger-dev/app
        cnci_tool = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/lnc_rna/CNCI/CNCI.py')
        self.logger.debug('soft {} cnci {}'.format(self.config.SOFTWARE_DIR, cnci_tool))
        cnci_out_dir = os.path.join(self.work_dir, 'cnci_workdir')

        cmd_fac_obj = CommandFactory(self.python_path + ' ' + cnci_tool)
        cmd_fac_obj.add_param('-f', self.option('fasta_file'), param_desc='fasta file path')
        cmd_fac_obj.add_param('-m', model, param_desc='classification models: plant --> pl, others --> ve')
        cmd_fac_obj.add_param('-o', cnci_out_dir, param_desc='output file path')
        cmd_fac_obj.add_param('-p', self.option('cpu'), param_desc='output file path')

        self.cmd_runner('cnci_predict', cmd_fac_obj.to_string(), check_stat=True)

        return os.path.join(cnci_out_dir, 'CNCI.index')

    def lncrna_mark(self, predict_file):
        # Transcript ID   index   score   start   end     length
        out_fields = ['transcript_id', 'transcript_length', 'orf_length', 'score', 'label']
        # extr_fields = ['Transcript ID', 'index', 'score', 'start', 'end', 'length']
        out_demo = '\t'.join('{%s}' % i for i in out_fields) + '\n'
        with open(predict_file) as in_handler, open(self.cnci_out_file, 'w') as out_handler:
            out_handler.write('\t'.join(out_fields) + '\n')
            for line_dic in csv.DictReader(in_handler, delimiter='\t'):
                t_id, gene_info = line_dic['Transcript ID'].split(' ')
                score = float(line_dic['score'])
                length = line_dic['length']
                try:
                    line_dic['orf_length'] = int(line_dic['end']) - int(line_dic['start'])
                except ValueError:
                    line_dic['orf_length'] = int(line_dic['length'])
                line_dic['label'] = 'noncoding' if score < self.option('cnci_score') else 'coding'
                out_handler.write(out_demo.format(
                    transcript_id=t_id, transcript_length=length, **line_dic))

    def run(self):
        super(CnciPredictorTool, self).run()
        predict_file = self.cnci_predict()
        self.lncrna_mark(predict_file)
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
                "id": "CnciPredictor_" + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "tool_lab.coding_prediction.cnci_predictor",
                "instant": False,
                "options": dict(
                    taxonmy="animal",
                    fasta_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab_batch2/Coding_predict/basic_filter.fa",
                    cnci_score=0,
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()

    unittest.main()
