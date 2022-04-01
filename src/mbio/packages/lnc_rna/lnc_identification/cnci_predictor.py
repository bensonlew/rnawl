#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/2/19 11:30
@file    : cnci_predictor.py
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
        self.cnci_score = 0
        self.cnci = r'/mnt/ilustre/users/isanger/app/bioinfo/lnc_rna/CNCI/CNCI.py'
        self.fa_file = None
        self.out_dir = './'
        self.out_file = 'cnci_output.txt'
        self.gtf_file = None
        self.process_num = 2
        self.is_ready = False
        self.__parser = argparse.ArgumentParser()
        self.__init_params()
        self.__parse_args()

    def __init_params(self):
        self.__parser.add_argument('-p', '--python_path', default='python', type=str)
        self.__parser.add_argument('-c', '--cnci_soft', default=self.cnci, type=str,
                                   help='program path, [default: %(default)s]')
        self.__parser.add_argument('-f', '--fasta_file', type=str,
                                   help='fasta file path, the seq head of fasta: ">transcrip_id"')
        self.__parser.add_argument('-g', '--gtf_file', type=str,
                                   help='gtf file which will be predicted')
        self.__parser.add_argument('-m', '--model', default=self.model, type=str, choices=['ve', 'pl'],
                                   help='classification models ("ve" for vertebrate species, '
                                        '"pl" for plat species), [default: %(default)s]')
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
        self.__parser.add_argument('-s', '--cnci_score', default=self.cnci_score, type=float,
                                   help='if score < cnci_score, it will be classified as noncoding, '
                                        '[default: %(default)s]')
        self.__parser.add_argument('-t', '--parallel', default=self.process_num, type=int,
                                   help='the number of parallel processes, [default: %(default)s]')

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
        self.cnci = self.__check_exist(args_obj.cnci_soft)
        self.fa_file = self.__check_exist(args_obj.fasta_file)
        self.gtf_file = self.__check_exist(args_obj.gtf_file)
        out_dir, self.out_file = os.path.split(args_obj.out_file)
        self.out_dir = self.__check_dir(out_dir)
        self.exon_number = args_obj.exon_number
        self.t_length = args_obj.trans_length
        self.orf_length = args_obj.orf_length
        self.cnci_score = args_obj.cnci_score
        self.process_num = args_obj.parallel
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


class CNCIPredictor(object):
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
                        'class_code': class_code
                    }
                    continue
                res_dic[transcript_id]['exon_number'] += 1
        sub_pipe.append(res_dic)

    def _stat_exon(self):
        pipe = multiprocessing.Manager().list()
        process = multiprocessing.Process(target=self._process_run, args=(pipe, self.param_obj.gtf_file))
        process.start()
        return process, pipe

    def _predict(self):
        temp_out_dir = os.path.join(self.param_obj.out_dir, 'temp_cnci')
        if not os.path.isdir(temp_out_dir):
            os.makedirs(temp_out_dir)
        else:
            os.system('rm %s/*' % temp_out_dir)
        cmd = '{python} {soft} -f {fa_file} -o {outdir} -p {proc_num} -m {model}'.format(
            python=self.param_obj.python,
            soft=self.param_obj.cnci,
            fa_file=self.param_obj.fa_file,
            outdir=temp_out_dir,
            proc_num=self.param_obj.process_num,
            model=self.param_obj.model,
        )

        process = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, universal_newlines=True)
        error_list = []
        for line in process.stderr:
            error_list.append(line)
        process.wait()
        if process.returncode != 0:
            raise Exception('%s -- ERROR: %s' % (cmd, ''.join(error_list)))
        out_file = os.path.join(temp_out_dir, 'CNCI.index')
        return out_file

    def _output(self, temp_out, add_dic):
        def res_filter(dic):
            class_codes = {'i', 'j', 'x', 'u', 'o'}
            coding_prob = dic['score']
            exon_num = dic['exon_number']
            class_code = dic['class_code']
            transcript_len = dic['length']
            orf_length = dic['ORF_length']
            if class_code in class_codes \
                    and float(coding_prob) < self.param_obj.cnci_score \
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
            'exon_number', 'class_code', 'cnci_score', 'cnci_label'
        ]
        # Transcript ID   index   score   start   end     length
        demo_list = [
            'transcript_id', 'gene_id', 'length',  'strand',
            'exon_number', 'class_code', 'score', 'index'
        ]
        demo_line = '\t'.join('{%s}' % i for i in demo_list) + '\n'
        out_file = os.path.join(self.param_obj.out_dir, self.param_obj.out_file)
        with open(temp_out) as infile, open(out_file, 'w') as outfile:
            outfile.write('\t'.join(head_list) + '\n')
            for line_dic in csv.DictReader(infile, delimiter='\t'):
                line_dic = {k.replace(' ', '_'): v for k, v in line_dic.items()}
                t_id = line_dic['Transcript_ID'].split(' ', 1)[0]
                t_dic = add_dic.get(t_id)
                line_dic['ORF_length'] = int(line_dic['end']) - int(line_dic['start'])
                # python3 format
                # outfile.write(demo_line.format(transcript_id=t_id, **line_dic, **t_dic))
                # python2 format
                line_dic.update(t_dic)
                line_dic['index'] = 'long_noncoding' if res_filter(line_dic) else line_dic['index']
                outfile.write(demo_line.format(transcript_id=t_id, **line_dic))

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
        cpc_obj = CNCIPredictor(param_obj=obj)
        cpc_obj.run()
    print(time.time() - s, ' S')
