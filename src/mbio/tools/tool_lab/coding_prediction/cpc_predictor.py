#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/2/26 12:51
@file    : cpc_predictor.py
"""
import csv
import os
import unittest

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
from mbio.packages.lnc_rna.lnc_identification.predict_common_tools import CommandFactory


class CpcPredictorAgent(Agent):
    def __init__(self, parent):
        super(CpcPredictorAgent, self).__init__(parent)
        # 'fasta_file', 'gtf_file', 'hexamer_dat', 'logit_model',
        # , 'transcript_len', 'cnci_score',
        # 'orf_len', 'cpat_score', 'taxonmy', 'exon_num'
        options = [
            {'name': 'fasta_file', 'type': 'string'},
            {'name': 'cpc_score', 'type': 'float', 'default': 0.5},
            {'name': 'cpu', 'type': 'int', 'default': 2}
        ]
        self.__cpu = 3
        self.add_option(options)
        self.option('cpu', self.__cpu)
        self.step.add_steps("cpc_predict")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.cpc_predict.start()
        self.step.update()

    def step_end(self):
        self.step.cpc_predict.finish()
        self.step.update()

    def check_options(self):
        fa_file = self.option('fasta_file')
        if fa_file and not os.path.isfile(fa_file):
            raise OptionError('缺少fasta file文件')
        return True

    def set_resource(self):
        self._cpu = 3
        self._memory = '5G'

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        super(CpcPredictorAgent, self).end()


class CpcPredictorTool(Tool):
    def __init__(self, config):
        super(CpcPredictorTool, self).__init__(config)
        self.python_path = "program/Python/bin/python"
        self.env_path = self.config.SOFTWARE_DIR + '/program/Python/bin'

        self.set_environ(PATH=self.env_path)
        self.cpc_out_file = os.path.join(self.output_dir, 'cpc_output.txt')

    def cmd_runner(self, cmd_name, cmd, check_stat=True, shell=False):
        cmd_obj = self.add_command(cmd_name, cmd, shell=shell)
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

    def cpc2_predict(self):
        # software_dir = /mnt/ilustre/users/sanger-dev/app
        cpc2_tool = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/lnc_rna/CPC2-beta/bin/CPC2.py')
        cpc_out_file = os.path.join(self.work_dir, 'temp_' + os.path.basename(self.cpc_out_file))
        self.logger.debug('soft {} cpc {}'.format(self.config.SOFTWARE_DIR, cpc2_tool))
        cmd_fac_obj = CommandFactory(self.python_path + ' ' + cpc2_tool)
        cmd_fac_obj.add_param('-i', self.option('fasta_file'))
        cmd_fac_obj.add_param('-o', cpc_out_file)
        self.cmd_runner('cpc2_predict', cmd_fac_obj.to_string(), check_stat=True)

        return cpc_out_file

    def lncrna_mark(self, cpc_file):
        # #ID  transcript_length  peptide_length  Fickett_score  pI  ORF_integrity  coding_probability  label
        out_fields = ['transcript_id', 'transcript_length', 'peptide_length',
                      'fickett_score', 'pi', 'orf_integrity', 'score', 'label']
        out_demo = '\t'.join('{%s}' % i for i in out_fields) + '\n'
        with open(cpc_file) as in_handler, open(self.cpc_out_file, 'w') as out_handler:
            out_handler.write('\t'.join(out_fields) + '\n')
            for line_dic in csv.DictReader(in_handler, delimiter='\t'):
                line_dic = {k.lower(): v for k, v in line_dic.items()}
                t_id = line_dic['#id']
                score = float(line_dic['coding_probability'])
                line_dic['label'] = 'noncoding' if score < self.option('cpc_score') else 'coding'
                out_handler.write(out_demo.format(transcript_id=t_id, score=score, **line_dic))

    def run(self):
        super(CpcPredictorTool, self).run()
        temp_cpc_file = self.cpc2_predict()
        self.lncrna_mark(temp_cpc_file)
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
                "id": "CpcPredictor_" + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "tool_lab.coding_prediction.cpc_predictor",
                "instant": False,
                "options": dict(
                    fasta_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab_batch2/Coding_predict/basic_filter.fa",
                    cpc_score=0.5,
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()

    unittest.main()
