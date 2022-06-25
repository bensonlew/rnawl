#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/2/26 16:44
@file    : pfam_scan_predictor.py
"""

import os
import re
import time
import unittest

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool

from mbio.files.lnc_rna.lnc_fasta import LncFastaFile
from mbio.packages.lnc_rna.lnc_identification.predict_common_tools import CommandFactory


class PfamScanPredictorAgent(Agent):
    def __init__(self, parent):
        super(PfamScanPredictorAgent, self).__init__(parent)
        # 'fasta_file', 'gtf_file', 'hexamer_dat', 'logit_model', 
        # , 'transcript_len', 'cnci_score',
        # 'orf_len', 'cpat_score', 'taxonmy', 'exon_num'
        options = [
            {'name': 'fasta_file', 'type': 'string', 'required': True},
            {'name': 'cpu', 'type': 'int', 'default': 12}
        ]
        self.add_option(options)
        self.step.add_steps("pfam_predict")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.pfam_predict.start()
        self.step.update()

    def step_end(self):
        self.step.pfam_predict.finish()
        self.step.update()

    def check_options(self):
        fa_file = self.option('fasta_file')
        if fa_file and not os.path.isfile(fa_file):
            raise OptionError('缺少fasta file文件')
        return True

    def set_resource(self):
        self._cpu = self.option('cpu')
        self.logger.debug('cpu: %s' % self._cpu)
        self._memory = '40G'

    def end(self):
        # self.add_upload_dir(self.output_dir)
        super(PfamScanPredictorAgent, self).end()


class PfamScanPredictorTool(Tool):
    def __init__(self, config):
        super(PfamScanPredictorTool, self).__init__(config)
        self.perl_rel_path = 'program/perl-5.24.0/bin/perl'
        env_path = ':'.join([
            os.path.join(self.config.SOFTWARE_DIR, 'miniconda2/bin'),
            os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/lnc_rna/hmmer/bin'),
            os.path.join(self.config.SOFTWARE_DIR, 'program/perl-5.24.0/bin')
        ])
        self.cmd_num = 1
        # export PERL5LIB, PERLLIB [/bioinfo/lnc_rna/PfamScan]
        self.pfam_tool_dir = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/lnc_rna/PfamScan')
        self.set_environ(PATH=env_path)
        self.set_environ(PERL5LIB=self.pfam_tool_dir, PERLLIB=self.pfam_tool_dir)
        self.pfam_work_dir = os.path.join(self.work_dir, 'pfam_work_dir')
        if not os.path.exists(self.pfam_work_dir):
            os.makedirs(self.pfam_work_dir)
        self.pfam_out_file = os.path.join(self.output_dir, 'pfam_output.txt')

    def cmd_runner(self, cmd_name, cmd, check_stat=True, shell=False):
        cmd_obj = self.add_command(cmd_name, cmd, shell=shell)
        if shell is True:
            cmd_obj.software_dir = ''
            cmd_obj._start_run_time = int(time.time())
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

    def obtain_all_rna(self):
        fasta_file = self.option('fasta_file')
        fa_heads = os.path.join(self.pfam_work_dir, 'fa_heads.txt')
        cmd_fac_obj = CommandFactory('/bin/grep')
        cmd_fac_obj.add_param('">"', fasta_file)
        cmd_fac_obj.add_param('>', fa_heads)

        self.cmd_runner('obtain_transcripts_list', cmd_fac_obj.to_string(), check_stat=True, shell=True)

        with open(fa_heads) as in_handler:
            ids_list = {item[1:].strip().split(' ')[0] for item in in_handler.read().strip().split('\n')}
        return ids_list

    def seq2protein(self):
        seq2protein = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/lnc_identification/seq2protein.pl')
        protein_out_file = os.path.join(self.pfam_work_dir, 'protein_output.txt')
        perl_path = os.path.join(self.config.SOFTWARE_DIR, self.perl_rel_path)
        cmd_fac_obj = CommandFactory(perl_path + ' ' + seq2protein)
        cmd_fac_obj.add_param('', self.option('fasta_file'), param_desc='fasta file path')
        cmd_fac_obj.add_param('>', protein_out_file, param_desc='out protein path')

        self.cmd_runner('pfam_seq2protein', cmd_fac_obj.to_string(), check_stat=True, shell=True)

        return protein_out_file

    def pfam_predict(self, pretein_fa, check_stat=True, pfam_outfile=None, batch_num=None):
        pfam_db = os.path.join(self.pfam_tool_dir, 'pfam_db')
        pfam_tool = os.path.join(self.pfam_tool_dir, 'pfam_scan.pl')
        if pfam_outfile is None:
            pfam_predict_file = os.path.join(self.pfam_work_dir, str(self.cmd_num) + '_pfam_predict_result.txt')
        else:
            pfam_predict_file = pfam_outfile

        # perl pfam_scan.pl -dir /mnt/ilustre/users/caiping.shi/software/PfamScan/pfam_db/
        # -outfile pfamscan_result -pfamB -fasta filter3.aa.fa -cpu 8
        cmd_fac_obj = CommandFactory(self.perl_rel_path + ' ' + pfam_tool)
        cmd_fac_obj.add_param('-dir', pfam_db, param_desc='pfam database')
        cmd_fac_obj.add_param('-outfile', pfam_predict_file, param_desc='out file')
        cmd_fac_obj.add_param('-fasta', pretein_fa, param_desc='pretein file')
        cmd_fac_obj.add_param('-cpu', '1', param_desc='cpu number')

        cmd_name = 'pfam_predict_' + str(self.cmd_num if batch_num is None else batch_num)
        cmd_obj = self.cmd_runner(cmd_name, cmd_fac_obj.to_string(), check_stat=check_stat)
        self.cmd_num += 1

        return pfam_predict_file, cmd_obj

    def lncrna_mark(self, *pfam_predict_files, **kwargs):
        handlers = (open(file) for file in pfam_predict_files)
        fields = None
        all_trancripts = set()
        signicant_transcripts = dict()
        ids_set = kwargs.get('ids_set', None)
        pattern = re.compile(' +')
        for in_handler in handlers:
            for line in in_handler:
                line = line.strip()
                if not line:
                    continue
                elif line.startswith('#'):
                    if fields or '>' not in line:
                        continue
                    fields = [item.strip().replace(' ', '_') for item in line.strip('<> \n').split('> <')]
                line_list = pattern.split(line.strip())
                t_id = line_list[0].rsplit('-', 1)[0]
                signicant = line_list[-2].strip()
                if signicant == '1':
                    signicant_transcripts[t_id] = 1
                if ids_set:
                    all_trancripts.add(t_id)
            # close file handler
            in_handler.close()

        all_trancripts = all_trancripts if ids_set is None else ids_set
        # output final result
        with open(self.pfam_out_file, 'w') as out_handler:
            out_handler.write('transcript_id\tlabel\n')
            out_list = []
            flag = 0
            for t in all_trancripts:
                out_list.append(t + ('\tcoding\n' if t in signicant_transcripts else '\tnoncoding\n'))
                flag += 1
                if flag % 1000 == 0:
                    out_handler.write(''.join(out_list))
                    out_list = []
            if out_list:
                out_handler.write(''.join(out_list))

    def split_faa(self, faa, total_num):
        def split_(faa_, split_sep_):
            flag = 0
            file_flag = 1
            out_handler = None
            out_files = []
            with LncFastaFile(file_path_=faa_) as fa_handler:
                for head, seq in fa_handler:
                    if flag >= split_sep_:
                        flag = 0
                    if flag == 0:
                        if out_handler is not None:
                            out_handler.close()
                        path = os.path.join(self.pfam_work_dir, 'split_' + str(file_flag) + '_proteins.faa')
                        out_files.append(path)
                        out_handler = open(path, 'w')
                        file_flag += 1
                    flag += 1
                    out_handler.write(">" + str(head) + "\n" + str(seq) + "\n")
            return out_files

        split_sep = 360
        is_split = False
        if total_num / split_sep > 10:
            split_sep = total_num / 10 + 1
            is_split = True
        elif total_num / split_sep > 1.5:
            split_sep = 360
            is_split = True
        if is_split:
            files = split_(faa_=faa, split_sep_=split_sep)
        else:
            files = [faa]

        return files

    def batch_run(self, *files):
        """wait until compcompletion

        :param files:
        :return:
        """
        run_num = self.option('cpu') - 1
        run_num = 1 if run_num == 0 else run_num
        predict_results = []
        start = 0
        end = run_num
        batch_num = 0
        while True:
            if start > len(files):
                break

            predict_results = []
            cmd_objs = []
            for file in files[start: end]:
                self.logger.debug(file)
                batch_num += 1
                pfam_outfile = os.path.join(self.pfam_work_dir, 'batch-%s_pfam_predict_result.txt' % batch_num)
                predict_file, cmd_obj = self.pfam_predict(file, check_stat=False, pfam_outfile=pfam_outfile,
                                                          batch_num=batch_num)
                predict_results.append(predict_file)
                cmd_objs.append(cmd_obj)
            self.check_stat(*cmd_objs)

            start = end
            end += run_num

        return predict_results

    def run(self):
        super(PfamScanPredictorTool, self).run()
        all_trans_ids = self.obtain_all_rna()
        # fa_seq to pretein seq
        pretein_fa = self.seq2protein()

        # 判断用于切分文件
        files = self.split_faa(faa=pretein_fa, total_num=len(all_trans_ids) * 3)
        # 批量预测
        predict_results = self.batch_run(*files)

        # label seq with coding and noncoding
        self.lncrna_mark(*predict_results, ids_set=all_trans_ids)
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
                "id": "PfamScanPredictor_" + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "tool_lab.coding_prediction.pfam_scan_predictor",
                "instant": False,
                "options": dict(
                    fasta_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab_batch2/Coding_predict/basic_filter.fa",
                    cpu=12
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()


    unittest.main()
