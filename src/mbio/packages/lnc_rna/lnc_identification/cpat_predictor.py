#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/2/19 15:08
@file    : cpat_predictor.py
@author  : zhipeng.zhao
"""

import argparse
import csv
import sys
import os
import subprocess
import multiprocessing


class ParamsObject(object):
    def __init__(self):
        self.model = 've'
        self.exon_number = 2
        self.t_length = 200
        self.orf_length = 300
        self.cpat_score = 0.5
        self.cpat = r'cpat.py'
        self.fa_file = None
        self.out_dir = './'
        self.out_file = 'cpat_output.txt'
        self.logit_model = None
        self.hexamer_dat = None
        self.gtf_file = None
        self.is_ready = False
        self.__parser = argparse.ArgumentParser()
        self.__init_params()
        self.__parse_args()

    def __init_params(self):
        self.__parser.add_argument('-p', '--python_path', default='python', type=str)
        self.__parser.add_argument('-c', '--cpat_soft', default=self.cpat, type=str,
                                   help='program path, [default: %(default)s]')
        self.__parser.add_argument('-f', '--fasta_file', type=str,
                                   help='fasta file path, the seq head of fasta: ">transcrip_id"')
        self.__parser.add_argument('-g', '--gtf_file', type=str,
                                   help='gtf file which will be predicted')
        self.__parser.add_argument('-m', '--logit_model', default=self.model, type=str, required=True,
                                   help='logit model out of your own training datset, '
                                        'Examples: .../dat/Human_logitModel.RData')
        self.__parser.add_argument('-x', '--hexamer_dat', type=str,
                                   help='this table out of your own training dataset, '
                                        'from running `make_hexamer_tab.py`'
                                        'Examples: .../dat/Human_Hexamer.tsv')
        self.__parser.add_argument('-o', '--out_file', default=self.out_file, type=str,
                                   help='output file path, [default: %(default)s]')
        self.__parser.add_argument('-e', '--exon_number', default=self.exon_number, type=int,
                                   help='the number of exon number in transcript used filter the '
                                        'needed transcript which exon number greater be equal or '
                                        'greater than `exon_number` [default: %(default)s]')
        self.__parser.add_argument('-l', '--trans_length', default=self.t_length, type=int,
                                   help='the minimum length of transcript length, '
                                        '[default: %(default)s]')
        self.__parser.add_argument('-r', '--orf_length', default=self.orf_length, type=int,
                                   help='the orf length of transcript, [default: %(default)s]')
        self.__parser.add_argument('-s', '--cpat_score', default=self.cpat_score, type=float,
                                   help='if score < cpat_score, it will be classified as noncoding, '
                                        '[default: %(default)s]')

    def __check_exist(self, path):
        if path is None or (path is not None and not os.path.isfile(path)):
            raise Exception('%s is not exist' % path)
        return path

    def __check_dir(self, path):
        if path is None:
            raise Exception('%s is not exist')
        elif path == '.' or path == '':
            return path
        elif not os.path.isdir(path):
            os.makedirs(path)
        return path

    def __parse_args(self):
        argv = sys.argv
        if len(argv) < 2 or '-h' in argv or '--help' in argv:
            self.__parser.parse_args(['-h'])

        args_obj = self.__parser.parse_args()
        self.python = args_obj.python_path \
            if args_obj.python_path == 'python' else self.__check_exist(args_obj.python_path)
        self.cpat = self.cpat if args_obj.cpat_soft == r'cpat.py' \
            else self.__check_exist(args_obj.cpat_soft)
        self.fa_file = self.__check_exist(args_obj.fasta_file)
        self.gtf_file = self.__check_exist(args_obj.gtf_file)
        out_dir, self.out_file = os.path.split(args_obj.out_file)
        self.out_dir = self.__check_dir(out_dir)
        self.logit_model = self.__check_exist(args_obj.logit_model)
        self.hexamer_dat = self.__check_exist(args_obj.hexamer_dat)
        self.exon_number = args_obj.exon_number
        self.t_length = args_obj.trans_length
        self.orf_length = args_obj.orf_length
        self.cpat_score = args_obj.cpat_score
        self.is_ready = True


class GTFParser(object):
    def __init__(self, file):
        self.__file = file
        self.__file_handler = None

    def __iter__(self):
        if self.__file_handler is None:
            self.__file_handler = open(self.__file)
        for line in self.__file_handler:
            line_splits  = line.strip().split('\t')
            line_splits[3] = int(line_splits[3])
            line_splits[4] = int(line_splits[4])
            line_splits[8] = {
                k.strip(): v.strip(' "')
                for k, v in (
                    item.strip().split(' ', 1) for item in line_splits[8].strip('; ').split(';')
                )
            }
            yield line_splits

    def __enter__(self):
        self.__file_handler = open(self.__file)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__file_handler.close()


class CPATPredictor(object):
    def __init__(self, param_obj):
        self.param_obj = param_obj

    @staticmethod
    def _process_run(sub_pipe, gtf_file):
        with GTFParser(gtf_file) as handler:
            res_dic = {}
            for line_items in handler:
                type_ = line_items[2]
                strand = line_items[6]
                prop_dic = line_items[8]
                gene_id = prop_dic['gene_id']
                transcript_id = prop_dic['transcript_id']
                if type_ == 'transcript':
                    res_dic.setdefault(transcript_id, {})
                    class_code = prop_dic['class_code']
                    res_dic[transcript_id] = {
                        'gene_id': gene_id,
                        'strand': strand,
                        'exon_number': 0,
                        'class_code': class_code,
                    }
                    continue
                start = int(line_items[3])
                end = int(line_items[4])
                res_dic[transcript_id]['exon_number'] += 1
        sub_pipe.append(res_dic)

    def _stat_exon(self):
        pipe = multiprocessing.Manager().list()
        process = multiprocessing.Process(target=self._process_run, args=(pipe, self.param_obj.gtf_file))
        process.start()
        return process, pipe

    def _predict(self):
        temp_file = os.path.join(self.param_obj.out_dir, 'temp_' + self.param_obj.out_file)
        # cpat.py -g temp.fa -d ~/app/bioinfo/lnc_rna/CPAT-1.2.4/dat/Human_logitModel.RData
        #       -x ~/app/bioinfo/lnc_rna/CPAT-1.2.4/dat/Human_Hexamer.tsv -o output2

        # cmd = '{python} {soft} -g {fa_file} -d {logit_model} -x {hexamer_dat} -o {outfile}'.format(
        cmd = '{soft} -g {fa_file} -d {logit_model} -x {hexamer_dat} -o {outfile}'.format(
            soft=self.param_obj.cpat,
            fa_file=self.param_obj.fa_file,
            logit_model=self.param_obj.logit_model,
            hexamer_dat=self.param_obj.hexamer_dat,
            outfile=temp_file
        )
        process = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, universal_newlines=True)
        error_list = []
        for line in process.stderr:
            error_list.append(line)
        process.wait()
        if process.returncode != 0:
            raise Exception('%s -- ERROR: %s' % (cmd, ''.join(error_list)))
        return temp_file

    def _output(self, temp_out, add_dic):
        def res_filter(dic):
            class_codes = {'i', 'j', 'x', 'u', 'o'}
            coding_prob = dic['coding_prob']
            exon_num = dic['exon_number']
            class_code = dic['class_code']
            transcript_len = dic['mRNA_size']
            orf_length = dic['ORF_size']
            if class_code in class_codes \
                    and float(coding_prob) < self.param_obj.cpat_score \
                    and int(exon_num) >= self.param_obj.exon_number \
                    and int(transcript_len) >= self.param_obj.t_length \
                    and int(orf_length) <= self.param_obj.orf_length:
                return True
            return False
        # head_list = [
        #     'transcript_id', 'gene_id', 'transcript_length', 'peptide_length',
        #     'ORF_length', 'strand', 'exon_number', 'Fickett_score', 'class_code',
        #     'pI', 'ORF_integrity', 'coding_probability', 'label'
        # ]
        head_list = [
            'transcript_id', 'gene_id', 'transcript_length',  'strand',
            'exon_number', 'class_code', 'cpat_score', 'cpat_label'
        ]
        # mRNA_size  ORF_size  Fickett_score  Hexamer_score  coding_prob
        demo_list = [
            'transcript_id', 'gene_id', 'mRNA_size',  'strand',
            'exon_number', 'class_code', 'coding_prob', 'label'
        ]
        demo_line = '\t'.join('{%s}' % i for i in demo_list) + '\n'
        out_file = os.path.join(self.param_obj.out_dir, self.param_obj.out_file)
        with open(temp_out) as infile, open(out_file, 'w') as outfile:
            outfile.write('\t'.join(head_list) + '\n')
            head_line = infile.readline().strip().split('\t')
            head_line.insert(0, 'transcript_id')
            for line in infile:
                l_split = line.strip().split('\t')
                line_dic = {k: v for k, v in zip(head_line, l_split)}
                line_dic = {k.replace(' ', '_'): v for k, v in line_dic.items()}
                t_id = line_dic['transcript_id'].split(' ', 1)[0]
                t_dic = add_dic.get(t_id)
                # python3 format
                # outfile.write(demo_line.format(transcript_id=t_id, **line_dic, **t_dic))
                # python2 format
                line_dic.update(t_dic)
                line_dic['label'] = ('long_noncoding' if res_filter(line_dic) else 'noncoding') \
                    if float(line_dic['coding_prob']) < self.param_obj.cpat_score else 'coding'
                outfile.write(demo_line.format(**line_dic))

    def run(self):
        process, exon_stat_list = self._stat_exon()
        temp_out = self._predict()
        process.join()
        self._output(temp_out, exon_stat_list[0])


if __name__ == '__main__':
    import time
    s = time.time()
    obj = ParamsObject()
    if obj.is_ready:
        cpc_obj = CPATPredictor(param_obj=obj)
        cpc_obj.run()
    print(time.time() - s, ' S')
