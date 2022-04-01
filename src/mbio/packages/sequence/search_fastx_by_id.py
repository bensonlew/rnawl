# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

import re


def search_fasta_by_id(fasta, fasta_id):
    try:
        with open(fasta, 'r') as f, open('searchID_result.fasta', 'w') as out_file:
            match = 0
            for line in f:
                for i in fasta_id.split(','):
                    # print i
                    if re.match(r'>', line) and i in line:
                        out_file.write('{}'.format(line))
                        out_file.write('{}'.format(f.next()))
                        match += 1
                        # print match
        return match
    except IOError:
        print '无法打开fasta文件'


def search_fastq_by_id(fastq, fastq_id):
    try:
        with open(fastq, 'r') as f, open('searchID_result.fastq', 'w') as out_file:
            match = 0
            for line in f:
                for i in fastq_id.split(','):
                    if re.match(r'@', line) and i in line:
                        out_file.write('{}'.format(line))
                        out_file.write('{}{}{}'.format(f.next(), f.next(), f.next()))
                        match += 1
        return match
    except IOError:
        print '无法打开fastq文件'


def search_fasta_by_idfile(fasta, id_file):
    try:
        with open(fasta, 'rb') as f, open(id_file, 'rb') as idfile, open('searchID_result.fasta', 'wb') as out_file:
            match = 0
            id_list = idfile.readlines()
            for line in f:
                for fas_id in id_list:
                    if re.match(r'>', line) and fas_id.strip() in line:
                        out_file.write('{}'.format(line))
                        out_file.write('{}'.format(f.next()))
                        match += 1
                        print match
        return match, id_list
    except IOError:
        print '无法打开fasta文件'


def search_fastq_by_idfile(fastq, id_file):
    try:
        with open(fastq, 'rb') as f, open(id_file, 'rb') as idfile, open('searchID_result.fastq', 'wb') as out_file:
            match = 0
            id_list = idfile.readlines()
            for line in f:
                for fas_id in id_list:
                    if re.match(r'@', line) and fas_id.strip() in line:
                        out_file.write('{}'.format(line))
                        out_file.write('{}{}{}'.format(f.next(), f.next(), f.next()))
                        match += 1
                        print match
        return match, id_list
    except IOError:
        print '无法打开文件'
