#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2019/3/15 12:56
@file    : ids_mapping.py
"""
import argparse
import itertools
import os
import sys
from collections import defaultdict


class FastaFile(object):
    def __init__(self, file):
        super(FastaFile, self).__init__()
        self.__file_handler = None
        self.seq_list = []
        self.seq_head = None
        self.file_path = file

    def __init_handler(self):
        self.__file_handler = open(self.file_path)
        while True:
            line = self.__file_handler.readline()
            if line == '\n' or line.startswith('#'):
                continue
            elif line.startswith('>'):
                self.seq_head = line
                break

    def __enter__(self):
        if self.__file_handler is None:
            self.__init_handler()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__file_handler.close()

    def __iter__(self):
        if self.__file_handler is None:
            self.__init_handler()
        return self

    def __next__(self):
        """
        python3 iterator iter method
        :return:
        """
        if self.__file_handler is None:
            self.__init_handler()
        while True:
            line = self.__file_handler.readline()
            if line.startswith('>') or (not line and len(self.seq_list) > 0):
                temp_list = self.seq_list
                temp_head = self.seq_head
                self.seq_list = []
                self.seq_head = line
                return temp_head, ''.join(temp_list)
            elif not line:
                raise StopIteration
            self.seq_list.append(line)

    def next(self):
        """
        python2 iterator iter method
        :return:
        """
        return self.__next__()


class Params(object):
    def __init__(self):
        self.fasta = None
        self.lnc_ids = None
        self.outdir = '.'
        self.ids_mapping = None
        self.ref_lnc_list = None
        self.ref_db = None
        self.db_name = None
        self.is_ready = False
        self.__args_check(self.__add_args())

    def __add_args(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('-f', '--fasta', type=str, help='fasta file path', required=True)
        parser.add_argument('-l', '--lnc_ids', type=str, help='lncRNA ids list')
        parser.add_argument('-o', '--outdir', type=str, help='output directory', default='.')
        ids_mapping = parser.add_argument_group('ids mapping arguments:',
                                                description='The following parameters must appear together')
        ids_mapping.add_argument('-i', '--ids_mapping', type=str, help='self id to ref id')
        ids_mapping.add_argument('-r', '--ref_db', type=str, help='ref lncRNA data database name')
        ids_mapping.add_argument('-d', '--db_name', type=str, help='database name of oneself')
        ids_mapping.add_argument('-n', '--ref_lnc_list', type=str, help='ref db lnc list to filter fa')

        return parser

    def __args_check(self, parser):
        argvs = sys.argv
        if len(argvs) < 3:
            parser.parse_args(['-h'])
        args_obj = parser.parse_args()
        # if all((args_obj.ids_mapping, args_obj.ref_db, args_obj.db_name, args_obj.ref_lnc_list)):
        #     self.ids_mapping = args_obj.ids_mapping
        #     self.ref_db = args_obj.ref_db
        #     self.db_name = args_obj.db_name
        #     self.ref_lnc_list = args_obj.ref_lnc_list
        if args_obj.ids_mapping:
            self.ids_mapping = args_obj.ids_mapping
        if args_obj.ref_db:
            self.ref_db = args_obj.ref_db
        if args_obj.db_name:
            self.db_name = args_obj.db_name
        if args_obj.ref_lnc_list:
            self.ref_lnc_list = args_obj.ref_lnc_list
        self.fasta = args_obj.fasta
        self.lnc_ids = args_obj.lnc_ids
        self.outdir = args_obj.outdir


        if not os.path.isdir(self.outdir):
            os.makedirs(self.outdir)

        self.is_ready = True
        return args_obj


class FastaFilter(object):

    def __init__(self, params_obj):
        self.params_obj = params_obj
        self.id_check = self._id_check()
        self.ref_lnc_map = self._ref_lnc_map()

    def _head_perser(self, head_line):
        t_id = head_line.strip().split()[0][1:]
        return t_id

    def _ref_lnc_map(self):
        lnc_set = set()
        if self.params_obj.ref_lnc_list and os.path.isfile(self.params_obj.ref_lnc_list):
            with open(self.params_obj.ref_lnc_list) as in_handler:
                lnc_set = {line.strip() for line in in_handler if line}
        def _filter(_id):
            return _id in lnc_set
        return _filter

    def _id_check(self):
        #ref_id_modify = (lambda i: i.strip().split('.')[0]) if self.params_obj.ref_db == 'ensembl' else (
        ref_id_modify = (lambda i: i.strip()) if self.params_obj.ref_db == 'ensembl' else (
            lambda line: line.strip())
        if self.params_obj.ids_mapping is not None:
            with open(self.params_obj.ids_mapping) as in_handler:
                ids_dict = {k.strip(): ref_id_modify(v) for k, v in (line.strip().split('\t') for line in in_handler)}
        else:
            ids_dict = {}

        if self.params_obj.lnc_ids is not None:
            with open(self.params_obj.lnc_ids) as in_handler:
                lnc_set = {line.strip() for line in in_handler}
                lnc_num = len(lnc_set)
        else:
            lnc_set = set()
            lnc_num = 0

        def __id_check(id_name):
            ref_id = None
            is_lncrna = False
            if lnc_num != 0 and id_name in lnc_set:
                ref_id = ids_dict.get(id_name)
                is_lncrna = True
            elif lnc_num == 0:
                ref_id = ids_dict.get(id_name)
                is_lncrna = True
            return is_lncrna, ref_id

        return __id_check

    def _run(self):
        basename = os.path.basename(self.params_obj.fasta)
        fa_out = os.path.join(self.params_obj.outdir, basename)
        map_ids = defaultdict(list)
        with FastaFile(self.params_obj.fasta) as fa_handler, open(fa_out, 'w') as out_handler:
            for head, seq in fa_handler:
                t_id = self._head_perser(head)
                if self.params_obj.db_name == 'ensembl':
                    # t_id = t_id.split('.')[0]
                    t_id = t_id
                is_lnc, ref_id = self.id_check(t_id)
                if is_lnc is True:
                    if ref_id is not None and self.ref_lnc_map(ref_id):
                        map_ids[ref_id].append(t_id)
                    elif ref_id is not None and self.ref_lnc_map(ref_id.split(".")[0]):
                        map_ids[ref_id.split(".")[0]].append(t_id)
                    else:
                        out_handler.write('>%s_%s\n' % (t_id, self.params_obj.db_name) + seq)

        total_len = len(map_ids)
        if total_len == 0:
            return

        mapped_ids_name = self.params_obj.db_name + '2' + self.params_obj.ref_db + '.txt'
        mapped_ids = os.path.join(self.params_obj.outdir, mapped_ids_name)
        with open(mapped_ids, 'w') as out_handler:
            head = 'transcript_id\t{}_transcript_id\n'.format(self.params_obj.db_name)
            out_handler.write(head)
            out_handler.write(''.join('{}\t{}\n'.format(k, ','.join(v)) for k, v in map_ids.items()))

    def run(self):
        self._run()


if __name__ == '__main__':
    params_obj = Params()
    if params_obj.is_ready:
        fa_obj = FastaFilter(params_obj)
        fa_obj.run()
