#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@time    : 2019/2/21 9:21
@file    : lncrna_predict.py
"""
import os
import unittest

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool


class LncrnaPredictAgent(Agent):
    def __init__(self, parent):
        super(LncrnaPredictAgent, self).__init__(parent)
        # 'fasta_file', 'gtf_file', 'hexamer_dat', 'logit_model', 
        # , 'transcript_len', 'cnci_score',
        # 'orf_len', 'cpat_score', 'taxonmy', 'exon_num'
        options = [
            {'name': 'taxonmy', 'type': 'string', 'default': 'Animal'},
            {'name': 'fasta_file', 'type': 'string'},
            {'name': 'gtf_file', 'type': 'string'},
            {'name': 'hexamer_dat', 'type': 'string', 'default': ''},
            {'name': 'logit_model', 'type': 'string', 'default': ''},
            {'name': 'transcript_len', 'type': 'int', 'default': 200},
            {'name': 'exon_num', 'type': 'int', 'default': 2},
            {'name': 'orf_len', 'type': 'int', 'default': 300},
            {'name': 'cpc_score', 'type': 'float', 'default': 0.5},
            {'name': 'cnci_score', 'type': 'float', 'default': 0},
            {'name': 'cpat_score', 'type': 'float', 'default': 0.5},
            {'name': 'identify_condition', 'type': 'int', 'default': 2},
            {'name': 'cpu', 'type': 'int', 'default': 2}
        ]
        cpu = 10
        self.add_option(options)
        self.option('cpu', cpu)
        self.step.add_steps("lncrna_predict")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.lncrna_predict.start()
        self.step.update()

    def step_end(self):
        self.step.lncrna_predict.finish()
        self.step.update()

    def check_options(self):
        fa_file = self.option('fasta_file')
        if fa_file or (fa_file and not os.path.isfile(fa_file)):
            OptionError('缺少fasta file文件')
        gtf_file = self.option('gtf_file')
        if gtf_file or (gtf_file and not os.path.isfile(gtf_file)):
            OptionError('缺少gtf file文件')

    def set_resource(self):
        self._cpu = 10
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(LncrnaPredictAgent, self).end()


class LncrnaPredictTool(Tool):
    def __init__(self, config):
        super(LncrnaPredictTool, self).__init__(config)
        self.python_path = "miniconda2/bin/python"
        self.perl_path = self.config.SOFTWARE_DIR + '/program/perl/perls/perl-5.24.0/bin/perl'
        gcc_path = self.config.SOFTWARE_DIR + '/gcc/8.1.0/bin'
        gcc_lib = self.config.SOFTWARE_DIR + '/gcc/8.1.0/lib:'
        gcc_lib += self.config.SOFTWARE_DIR + '/gcc/8.1.0/lib64'
        # /bioinfo/lnc_rna/PfamScan/pfam_scan.pl
        self.pfam_scan_dir = self.config.SOFTWARE_DIR + '/bioinfo/lnc_rna/PfamScan'
        # self.env_path = self.config.SOFTWARE_DIR + "/program/Python/bin:$PATH"
        # 不需要再bin后加:$PARH, 因为set_environ中自动添加了
        # 添加hmmer路径[ 目的使用hmmscan ]/bioinfo/lnc_rna/hmmer-3.2.1/bin
        self.env_path = self.config.SOFTWARE_DIR + '/program/Python/bin'
        self.env_path += self.config.SOFTWARE_DIR + '/bioinfo/lnc_rna/hmmer-3.2.1/bin'
        self.env_path += gcc_path
        # export PERL5LIB=/mnt/ilustre/users/isanger/app/bioinfo/lnc_rna/PfamScan:$PERL5LIB
        self.pfam_perl_lib = self.config.SOFTWARE_DIR + '/bioinfo/lnc_rna/PfamScan'
        self.set_environ(PATH=self.env_path, PERL5LIB=self.pfam_perl_lib, LD_LIBRARY_PATH=gcc_lib)
        self.cpc_out_file = None
        self.cnci_out_file = None
        self.cpat_out_file = None

    class CommandFactory(object):
        def __init__(self, exe_path):
            """
            :param exe_path: python cpat.py
            """
            self.__params_list = [exe_path]
            self.__cmd = None

        def add_param(self, param, param_value, param_desc=None):
            self.__params_list.extend((param, param_value))

        def to_string(self):
            if self.__cmd is None:
                self.__cmd = ' '.join(str(i) for i in self.__params_list)
            return self.__cmd

        def __str__(self):
            return self.to_string()

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

    def cpc2_predict(self):
        cpc2_pkg = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/cpc2_predictor.py')
        self.cpc_out_file = os.path.join(self.output_dir, 'cpc2_output.txt')

        cmd_fac_obj = self.CommandFactory(self.python_path + ' ' + cpc2_pkg)
        cmd_fac_obj.add_param('-f', self.option('fasta_file'))
        cmd_fac_obj.add_param('-g', self.option('gtf_file'))
        cmd_fac_obj.add_param('-o', self.cpc_out_file)
        cmd_fac_obj.add_param('-e', self.option('exon_num'), param_desc='minimum exon number in transcript')
        cmd_fac_obj.add_param('-l', self.option('transcript_len'), param_desc='minimum transcript length')
        cmd_fac_obj.add_param('-r', self.option('orf_len'), param_desc='maximum orf length')
        cmd_fac_obj.add_param('-s', self.option('cpc_score'), param_desc='< cnci_score --> lncRNA')

        return self.cmd_runner('cpc2_predict', cmd_fac_obj.to_string(), check_stat=False)

    def cnci_predict(self):
        taxonmy = self.option('taxonmy')
        model = 'pl' if taxonmy.lower() == 'plant' else 've'
        cnci_pkg = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/cnci_predictor.py')
        self.cnci_out_file = os.path.join(self.output_dir, 'cnci_output.txt')

        cmd_fac_obj = self.CommandFactory(self.python_path + ' ' + cnci_pkg)
        cmd_fac_obj.add_param('-f', self.option('fasta_file'), param_desc='fasta file path')
        cmd_fac_obj.add_param('-g', self.option('gtf_file'), param_desc='gtf file path')
        cmd_fac_obj.add_param('-m', model, param_desc='classification models: plant --> pl, others --> ve')
        cmd_fac_obj.add_param('-o', self.cnci_out_file, param_desc='output file path')
        cmd_fac_obj.add_param('-e', self.option('exon_num'), param_desc='minimum exon number in transcript')
        cmd_fac_obj.add_param('-l', self.option('transcript_len'), param_desc='minimum transcript length')
        cmd_fac_obj.add_param('-r', self.option('orf_len'), param_desc='maximum orf length')
        cmd_fac_obj.add_param('-s', self.option('cnci_score'), param_desc='< cnci_score --> lncRNA')
        cmd_fac_obj.add_param('-t', self.option('cpu'), param_desc='parallel num')

        return self.cmd_runner('cnci_predict', cmd_fac_obj.to_string(), check_stat=False)

    def cpat_predict(self):
        cnci_pkg = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/cpat_predictor.py')
        self.cpat_out_file = os.path.join(self.output_dir, 'cpat_output.txt')

        cmd_fac_obj = self.CommandFactory(self.python_path + ' ' + cnci_pkg)
        cmd_fac_obj.add_param('-f', self.option('fasta_file'), param_desc='fasta file path')
        cmd_fac_obj.add_param('-g', self.option('gtf_file'), param_desc='gtf file path')
        # logit model example: Human_logitModel.RData
        cmd_fac_obj.add_param('-m', self.option('logit_model'), param_desc='logit model')
        # hexamer dat, example: Human_Hexamer.tsv
        cmd_fac_obj.add_param('-x', self.option('hexamer_dat'), param_desc='hexamer dat')
        cmd_fac_obj.add_param('-o', self.cpat_out_file, param_desc='output file path')
        cmd_fac_obj.add_param('-e', self.option('exon_num'), param_desc='minimum exon number in transcript')
        cmd_fac_obj.add_param('-l', self.option('transcript_len'), param_desc='minimum transcript length')
        cmd_fac_obj.add_param('-r', self.option('orf_len'), param_desc='maximum orf length')
        cmd_fac_obj.add_param('-s', self.option('cpat_score'), param_desc='< cpat_score --> lncRNA')

        return self.cmd_runner('cpat_predict', cmd_fac_obj.to_string(), check_stat=False)

    def pfam_predict(self):
        pfam_scan_db = os.path.join(self.pfam_scan_dir, )

    def set_output(self):
        pass

    def run(self):
        super(LncrnaPredictTool, self).run()
        wait_list = []
        wait_list.append(self.cpc2_predict())
        wait_list.append(self.cnci_predict())
        if self.cpc_out_file is not None:
            wait_list.append(self.cpat_predict())
        self.check_stat(*wait_list)
        self.set_output()
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
                "id": "LncrnaPredict_" + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "lnc_rna.lncrna_predict",
                "instant": False,
                "options": dict(
                    taxonmy="animal",
                    fasta_file="/mnt/ilustre/users/sanger-dev/sg-users/zhaozhipeng/temp/test.fa",
                    gtf_file="/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/assemble/new_transcripts.gtf",
                    hexamer_dat="/mnt/ilustre/users/isanger/app/bioinfo/lnc_rna/CPAT-1.2.4/dat/Human_Hexamer.tsv",
                    logit_model="/mnt/ilustre/users/isanger/app/bioinfo/lnc_rna/CPAT-1.2.4/dat/Human_logitModel.RData",
                    transcript_len=200,
                    exon_num=2,
                    orf_len=300,
                    cpc_score=0.5,
                    cnci_score=0,
                    cpat_score=0.5,
                    identify_condition=2
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()

    unittest.main()
