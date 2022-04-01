#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/2/26 8:51
@file    : basic_filter.py
"""
import json
import os
import time
import unittest

import pandas as pd

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
from mbio.packages.lnc_rna.lnc_identification.predict_common_tools import CommandFactory, GTFParser, FastaParser


class BasicFilterAgent(Agent):
    def __init__(self, parent):
        super(BasicFilterAgent, self).__init__(parent)
        # 'fasta_file', 'gtf_file', 'hexamer_dat', 'logit_model',
        # , 'transcript_len', 'cnci_score',
        # 'orf_len', 'cpat_score', 'taxonmy', 'exon_num'
        options = [
            {'name': 'fasta_file', 'type': 'string'},
            {'name': 'gtf_file', 'type': 'string'},
            {'name': 'transcript_len', 'type': 'int', 'default': 200},
            {'name': 'exon_num', 'type': 'int', 'default': 2},
            {'name': 'orf_len', 'type': 'int', 'default': 300},
            {'name': 'out_file', 'type': 'string', 'default': 'basic_filter_1.fa'},
            {'name': 'known_transcripts_json', 'type': 'string'},
            {'name': 'cpu', 'type': 'int', 'default': 2}
        ]
        cpu = 10
        self.add_option(options)
        self.option('cpu', cpu)
        self.step.add_steps("basic_filter")
        self.on("start", self.step_start)
        self.on("end", self.step_end)

    def step_start(self):
        self.step.basic_filter.start()
        self.step.update()

    def step_end(self):
        self.step.basic_filter.finish()
        self.step.update()

    def check_options(self):
        fa_file = self.option('fasta_file')
        if not fa_file or (fa_file and not os.path.isfile(fa_file)):
            raise OptionError('缺少fasta file文件')
        gtf_file = self.option('gtf_file')
        if not gtf_file or (gtf_file and not os.path.isfile(gtf_file)):
            raise OptionError('缺少gtf file文件')
        json_file = self.option('known_transcripts_json')
        if not gtf_file or (json_file and not os.path.isfile(json_file)):
            raise OptionError('缺少 known_transcripts_json 文件')
        return True

    def set_resource(self):
        self._cpu = 10
        self._memory = '10G'

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        super(BasicFilterAgent, self).end()


class BasicFilterTool(Tool):
    def __init__(self, config):
        super(BasicFilterTool, self).__init__(config)
        env_path = ':'.join([self.config.SOFTWARE_DIR + '/program/Python/bin',
                            self.config.SOFTWARE_DIR + '/program/perl/perls/perl-5.24.0/bin/'])
        self.orf_stat_file = None
        self.set_environ(PATH=env_path)
        self.filter_work_dir = self.__dir_check(os.path.join(self.work_dir, 'basic_filter'))

    def __dir_check(self, dir_path):
        if not os.path.isdir(dir_path):
            os.makedirs(dir_path)
        return dir_path

    def __cmd_runner(self, cmd_name, cmd, check_stat=True, shell=False):
        cmd_obj = self.add_command(cmd_name, cmd, shell=shell)
        if shell is True:
            cmd_obj._start_run_time = int(time.time())
            cmd_obj.software_dir = ''
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

    def longorf_extract(self):
        perl_path = self.config.SOFTWARE_DIR + '/program/perl/perls/perl-5.24.0/bin/perl'
        tool_path = os.path.join(self.config.PACKAGE_DIR, 'lnc_rna/lnc_identification/extract_longorf.pl')
        self.orf_stat_file = os.path.join(self.output_dir, 'orf_stat.txt')
        cmd_fac_obj = CommandFactory(perl_path + ' ' + tool_path)
        cmd_fac_obj.add_param('--input', self.option('fasta_file'))
        cmd_fac_obj.add_param('>', self.orf_stat_file)

        self.__cmd_runner('extract_longorf', cmd_fac_obj.to_string(), shell=True, check_stat=True)

        df = pd.read_table(self.orf_stat_file, sep='\t', header=None, index_col=None)
        return {k: v for k, v in zip(df.ix[:, 0], df.ix[:, 1])}

    def gtf_filter(self, orf_dic):
        trans_set = set()
        class_codes = {'i', 'j', 'x', 'u', 'o'}
        def _filter(recored):
            transcript_len = recored['transcript_len']
            exon_num = recored['exon_num']
            class_code = recored['class_code']
            orf_len = recored['orf_len']
            if transcript_len >= self.option('transcript_len') \
                    and exon_num >= self.option('exon_num') \
                    and orf_len <= self.option('orf_len') \
                    and class_code in class_codes:
                trans_set.add(recored['id'])
            recored.clear()

        with GTFParser(self.option('gtf_file')) as gtf_handler:
            transcript_id = None
            one_record = {}
            t2g_ids_map = {}
            for line_list, line in gtf_handler:
                record_type = line_list[2]
                if record_type == 'exon':
                    attr_dic = line_list[8]
                    t_id = attr_dic.get('transcript_id')
                    if t_id not in t2g_ids_map:
                        t2g_ids_map[t_id] = {
                            'strand': line_list[6],
                            'gene_id': attr_dic.get('gene_id'),
                            'length': 0
                        }
                    if transcript_id and t_id != transcript_id:
                        _filter(one_record)
                    if len(one_record) == 0:
                        transcript_id = t_id
                        one_record['id'] = t_id
                        one_record['transcript_len'] = 0
                        one_record['exon_num'] = 0
                        one_record['class_code'] = attr_dic.get('class_code')
                        one_record['orf_len'] = orf_dic.get(t_id)
                    exon_len = line_list[4] - line_list[3] + 1
                    t2g_ids_map[t_id]['length'] += exon_len
                    one_record['transcript_len'] += exon_len
                    one_record['exon_num'] += 1
            _filter(one_record)
        self.out_ids_map(t2g_ids_map)
        return trans_set

    def fasta_filter(self, transcript_set, known_trans_set):
        self.filtered_fasta_file = os.path.join(self.output_dir, self.option('out_file'))
        raw_mrna_fasta_file = os.path.join(self.output_dir, 'raw_mrna.fa')
        head_parser = lambda line: {'id': line[1:].strip().split(' ')[0], 'head_line': line}
        with FastaParser(self.option('fasta_file'), record_head_parser=head_parser) as fa_handler, \
                open(raw_mrna_fasta_file, 'w') as mrna_out_handler, \
                open(self.filtered_fasta_file, 'w') as out_handler:
            for record_dic in fa_handler:
                t_id = record_dic['head']['id']
                if t_id in known_trans_set:
                    continue
                elif t_id in transcript_set:
                    out_handler.write(record_dic['head']['head_line'] + record_dic['fa_seq'])
                else:
                    mrna_out_handler.write(record_dic['head']['head_line'] + record_dic['fa_seq'])

    def known_trans_set(self):
        with open(self.option('known_transcripts_json')) as in_handler:
            return {k for k in json.load(in_handler)}

    def out_ids_map(self, ids_dict):
        """
        dump ids map to json file
        :return:
        """
        out_file = os.path.join(self.output_dir, 'transcript2gene_dict.json')
        with open(out_file, 'w') as out_handler:
            json.dump(ids_dict, out_handler, indent=4)

    # def set_output(self):
    #
    #     link = os.path.join(self.output_dir, os.path.basename(self.filtered_fasta_file))
    #     if os.path.exists(link):
    #         os.remove(link)
    #     os.link(self.filtered_fasta_file, link)

    def run(self):
        super(BasicFilterTool, self).run()
        orf_dic = self.longorf_extract()
        k_trans_set = self.known_trans_set()
        trans_set = self.gtf_filter(orf_dic)
        self.fasta_filter(trans_set, k_trans_set)
        # self.set_output()
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
                "id": "BasicFilter_" + str(random.randint(1, 10000)),
                "type": "tool",
                "name": "lnc_rna.lncrna_identification.basic_filter",
                "instant": False,
                "options": dict(
                    fasta_file="/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/assemble/"
                               "method-cufflinks/output/NewTranscripts/new_transcripts.fa",
                    gtf_file="/mnt/ilustre/users/isanger/sg-users/qinjincheng/lnc_rna/assemble/"
                             "method-cufflinks/output/NewTranscripts/new_transcripts.gtf",
                    known_transcripts_json=r"/mnt/ilustre/users/sanger-dev/workspace/20190410/"
                                           r"Single_extr_known_lnc_9325_7266/ExtrKnownLnc/output/known_lnc_in_new.json",
                    transcript_len=200,
                    exon_num=2,
                    orf_len=300,
                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()


    unittest.main()
