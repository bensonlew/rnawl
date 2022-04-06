## !/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/python
# -*- coding: utf-8 -*-
# __author__ = "yitong.feng"
# 20180719

import argparse
import os
import re
from collections import defaultdict
import shutil
from biocluster.config import Config


class rock_index(object):
    def __init__(self, input_file, input_type, fna):
        self.input_file = input_file
        self.input_type = input_type
        self.fna = fna
        self.chr_names = self.get_chr_names(self.input_file)

    def parse_each_chromosome_feature(self, in_file):
        f1 = open(in_file)
        header = f1.readline().strip().split('\t')
        f2 = open(in_file.split('.feature_table')[0] + '.ptt', 'w')
        f2.write('Infomation from: ' + in_file + '\n')
        f2.write('[int] proteins\n')
        f2.write('Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n')
        f3 = open(in_file.split('.feature_table')[0] + '.rnt', 'w')
        f3.write('Infomation from: ' + in_file + '\n')
        f3.write('[int] no_coding rnas\n')
        f3.write('Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n')
        for line in f1:
            tmp = line.strip('\n').split('\t')
            if not tmp[14]:
                tmp[14] = '-'
            if not tmp[10]:
                tmp[10] = '-'

            nt_len = int(tmp[8]) - int(tmp[7])
            pt_len = nt_len // 3
            if tmp[0] == 'CDS':
                if tmp[13].strip():
                    f2.write(tmp[7] + '..' + tmp[8] + '\t' + tmp[9] + '\t' + str(pt_len) +
                         '\t' + tmp[10] + '\t' + tmp[14] + '\t' + tmp[16] + '\t' +
                         '-' + '\t' + '-' + '\t' + tmp[13] + '\n')

            if (tmp[0] == 'rRNA') or (tmp[0] == 'tRNA'):
                if tmp[13].strip():
                    f3.write(tmp[7] + '..' + tmp[8] + '\t' + tmp[9] + '\t' + str(nt_len + 1) +
                         '\t' + tmp[10] + '\t' + tmp[14] + '\t' + tmp[16] + '\t' +
                         '-' + '\t' + '-' + '\t' + tmp[13] + '\n')
        f3.close()
        f2.close()
        f1.close()

    def parse_each_chromosome_gff(self, in_file):
        f1 = open(in_file)
        f2 = open(in_file.split('.gff')[0] + '.ptt', 'w')
        f2.write('Infomation from: ' + in_file + '\n')
        f2.write('[int] proteins\n')
        f2.write('Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n')
        f3 = open(in_file.split('.gff')[0] + '.rnt', 'w')
        f3.write('Infomation from: ' + in_file + '\n')
        f3.write('[int] no_coding rnas\n')
        f3.write('Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n')
        loca_tmp = ' '
        attr_tmp = ' '
        pseu_list = list()
        synonym_list = set()
        gfflines = f1.readlines()
        gene2product = dict()
        gene2pid = dict()
        rna_list = list()
        for line in gfflines:
            if not line.strip().startswith('#'):
                line = line.strip().split('\t')
                if len(line) >= 9:
                    name, feature, start, end, strand, attribute = line[0], line[2], line[3], line[4], line[6], line[8]
                    if feature.lower() == 'cds':
                        for att in attribute.split(';'):
                            if re.match('Parent=(.*)', att):
                                synonym = re.match('Parent=(.*)', att).group(1)
                            if re.match('^product=?"?(.*)"?$', att):
                                product = re.match('^product=?"?(.*)"?$', att).group(1)
                            if re.match('^Dbxref=.*\:(.*)$', att):
                                pid = re.match('^Dbxref=.*\:(.*)$', att).group(1)
                        try:
                            gene2product[synonym] = product
                            gene2pid[synonym] = pid
                        except:
                            pass
                    if feature == 'rRNA' or feature == 'tRNA':
                        for att in attribute.split(';'):
                            if re.match('gene_id "(.*)"', att):
                                synonym = re.match('gene_id "(.*)"', att).group(1)
                            if re.match('ID=(.*)', att):
                                synonym = re.match('ID=(.*)', att).group(1)
                                gene_ = re.match('ID=(.*)', att).group(1)
                            if re.match('locus_tag=(.*)', att):
                                synonym = re.match('locus_tag=(.*)', att).group(1)
                        if synonym:
                            rna_list.append(synonym)

        for line in gfflines:
            # synonym_list = set()
            name = feature = start = end = strand = attribute = ''
            gene_ = gene = pid = synonym = product = '-'
            if not line.strip().startswith('#'):
                line = line.strip().split('\t')
                if len(line) >= 9:
                    name, feature, start, end, strand, attribute = line[0], line[2], line[3], line[4], line[6], line[8]
                    location = line[3] + '..' + line[4]
                    ##----判断出pseudogene后会跟着俩重复基因的情况------
                    if feature == 'pseudogene':
                        if loca_tmp:
                            loca_tmp_ = loca_tmp
                            loca_tmp = location
                        else:
                            loca_tmp_ = loca_tmp = location
                        if loca_tmp_.startswith(line[3]) and not loca_tmp_.endswith(line[4]):
                            location = loca_tmp
                        if not loca_tmp_.startswith(line[3]) and loca_tmp_.endswith(line[4]):
                            continue
                    if feature == 'pseudogene':
                        if attr_tmp:
                            attr_tmp_ = attr_tmp
                            attr_tmp = attribute.split(';')[0]
                        else:
                            attr_tmp_ = attr_tmp = attribute.split(';')[0]
                        if attr_tmp_ == attribute.split(';')[0]:
                            continue
                        if attr_tmp in pseu_list:
                            continue
                        else:
                            pseu_list.append(attr_tmp)
                    # if loca_tmp.startswith(line[3]) and not loca_tmp.endswith(line[4]):
                    #     location = loca_tmp
                    # if not loca_tmp.startswith(line[3]) and loca_tmp.endswith(line[4]):
                    #     continue
                    for att in attribute.split(';'):
                        if re.match('gene_id "(.*)"', att):
                            synonym = re.match('gene_id "(.*)"', att).group(1)
                        if re.match('ID=(.*)', att):
                            synonym = re.match('ID=(.*)', att).group(1)
                            gene_ = re.match('ID=(.*)', att).group(1)
                        if re.match('locus_tag=(.*)', att):
                            synonym = re.match('locus_tag=(.*)', att).group(1)
                        if re.match('protein_id=?"?(.*)"?$', att):
                            pid = re.match('protein_id=?"?(.*)"?$', att).group(1)
                        if re.match('gene="?(.*)"?$', att):
                            gene = re.match('gene="?(.*)"?$', att).group(1)
                        if gene == "-" and re.match('Name="?(.*)"?$', att):
                            gene = re.match('Name="?(.*)"?$', att).group(1)
                        if gene == "-" and re.match('gene_name="?(.*)"?$', att):
                            gene = re.match('gene_name="?(.*)"?$', att).group(1)    # added by zhangyitong on 20210825
                    if gene_ in gene2product:
                        product = gene2product[gene_]
                        if not product:
                            product = '-'
                    if gene_ in gene2pid:
                        pid = gene2pid[gene_]
                    if synonym != '-':
                        if synonym in synonym_list and feature == 'CDS':
                            continue
                        else:
                            synonym_list.add(synonym)
                    # if feature == 'CDS':
                    if feature.lower() in ['gene', 'pseudogene'] and synonym not in rna_list and synonym != '-':
                        # f2.write(location + '\t' + strand + '\t' + str(int(end) - int(start) + 1) + '\t' + pid + '\t' + gene + '\t' + synonym + '\t' + '-\t' + '-\t' + product + '\n')
                        f2.write(location + '\t' + strand + '\t' + str((int(end) - int(start)) // 3) + '\t' + pid + '\t' + gene + '\t' + synonym + '\t' + '-\t' + '-\t' + product + '\n')
                    if feature == 'rRNA' or feature == 'tRNA' and synonym != '-':
                        # f3.write(location + '\t' + strand + '\t' + str(int(end) - int(start) + 1) + '\t' + pid + '\t' + gene + '\t' + synonym + '\t' + '-\t' + '-\t' + product + '\n')
                        f3.write(location + '\t' + strand + '\t' + str(int(end) - int(start)) + '\t' + pid + '\t' + gene + '\t' + synonym + '\t' + '-\t' + '-\t' + product + '\n')
        f3.close()
        f2.close()
        f1.close()

    def parse_each_chromosome_gtf(self, in_file):
        f1 = open(in_file)
        f2 = open(in_file.split('.gtf')[0] + '.ptt', 'w')
        f2.write('Infomation from: ' + in_file + '\n')
        f2.write('[int] proteins\n')
        f2.write('Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n')
        f3 = open(in_file.split('.gtf')[0] + '.rnt', 'w')
        f3.write('Infomation from: ' + in_file + '\n')
        f3.write('[int] no_coding rnas\n')
        f3.write('Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n')
        gtf_r = f1.readlines()
        cds_locations = list()
        for line in gtf_r:
            line = line.strip().split('\t')
            if len(line) >= 4:
                if line[2] == 'CDS':
                    cds_locations.append(line[3] + '..' + line[4])
        for line in gtf_r:
            line = line.strip().split('\t')
            name = feature = start = end = strand = attribute = ''
            gene = synonym = product = '-'
            if len(line) >= 8:
                name, feature, start, end, strand, attribute = line[0], line[2], line[3], line[4], line[6], line[8]
                location = line[3] + '..' + line[4]
                for att in attribute.split(';'):
                    if re.match('gene_name "(.*?)"', att.strip()):
                        gene = re.match('gene_name "(.*?)"', att.strip()).group(1)
                    if re.match('gene_id "(.*?)"', att.strip()):
                        synonym = re.match('gene_id "(.*?)"', att.strip()).group(1)
                if feature == 'CDS':
                    f2.write(location + '\t' + strand + '\t' + str(int(end) - int(start) + 1) + '\t' + '-' + '\t' + gene + '\t' + synonym + '\t' + '-\t' + '-\t' + product + '\n')
                if feature != 'CDS' and location not in cds_locations:
                    f3.write(location + '\t' + strand + '\t' + str(int(end) - int(start) + 1) + '\t' + '-' + '\t' + gene + '\t' + synonym + '\t' + '-\t' + '-\t' + product + '\n')
        f3.close()
        f2.close()
        f1.close()

    def get_chr_names(self, in_file):
        f1 = open(in_file)
        if self.input_type.lower() == 'feature':
            header = f1.readline()
            chr_names = set()
            for line in f1:
                tmp = line.strip().split('\t')
                chr_names.add(tmp[6])
        else:
            chr_names = set()
            for line in f1:
                if not line.startswith('#'):
                    tmp = line.strip().split('\t')
                    chr_names.add(tmp[0])
        f1.close()
        return chr_names

    def parse_feature_table(self, in_file):
        for chr_name in self.chr_names:
            f1 = open(in_file)
            f2 = open(chr_name + '.feature_table', 'w')
            f2.write(f1.readline())
            for line in f1:
                tmp = line.strip('\n').split('\t')
                if tmp[6] == chr_name:
                    f2.write(line)
            f2.close()
            f1.close()
            self.parse_each_chromosome_feature(chr_name + '.feature_table')

    def parse_gff(self, in_file):
        for chr_name in self.chr_names:
            f1 = open(in_file)
            f2 = open(chr_name + '.gff', 'w')
            for line in f1:
                tmp = line.strip().split('\t')
                if not tmp[0].startswith('#') and tmp[0] == chr_name:
                    f2.write(line)
            f2.close()
            f1.close()
            self.parse_each_chromosome_gff(chr_name + '.gff')

    def parse_gtf(self, in_file):
        for chr_name in self.chr_names:
            f1 = open(in_file)
            f2 = open(chr_name + '.gtf', 'w')
            for line in f1:
                tmp = line.strip('\n').split('\t')
                if not tmp[0].startswith('#') and tmp[0] == chr_name:
                    f2.write(line)
            f2.close()
            f1.close()
            self.parse_each_chromosome_gtf(chr_name + '.gtf')

    def split_fna(self, fna):
        chr_names = list()
        with open(fna) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                if line.startswith('>'):
                    chr_name = line.lstrip('>').split()[0]
                    chr_names.append(chr_name)
                    try:
                        f2.close()
                    except:
                        pass
                    finally:
                        f2 = open(chr_name + '.fna', 'w')
                        f2.write('>gi|xxx|ref|' + chr_name + '| ' + line[1:])
                    continue
                f2.write(line)
            f2.close()
        return chr_names

    def run(self):
        if self.input_type.lower() == 'feature':
            self.parse_feature_table(self.input_file)
        elif self.input_type.lower() == 'gff':
            self.parse_gff(self.input_file)
        elif self.input_type.lower() == 'gtf':
            self.parse_gtf(self.input_file)
        if self.fna:
            self.split_fna(self.fna)
        if os.path.exists('rock_index'):
            shutil.rmtree('rock_index')
        os.mkdir('rock_index')
        for chr_name in self.chr_names:
            os.mkdir('rock_index/' + chr_name)
            os.system('mv {}.* {}'.format(chr_name, 'rock_index/' + chr_name))


if __name__ == '__main__':
    """
    format ptt and rnt
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', required=True, help='feature table file or gff|gtf')
    parser.add_argument('-fna', default=None, help='genome fasta file to be parsed')
    parser.add_argument('-type', default='feature', help='can be feature,gff or gtf')
    args = parser.parse_args()
    # -----过滤gtf或gff超出参考基因组的情况-----
    if args.type != 'feature':
        bioawk_path = Config().SOFTWARE_DIR + '/bioinfo/seq/bioawk/'
        bedtools_path = Config().SOFTWARE_DIR + '/bioinfo/seq/bedtools-2.25.0/bin/'
        script_path = Config().SOFTWARE_DIR + '/bioinfo/rna/scripts/'
        cmd_bed = "bash " + script_path + "fasta_range.sh %s %s %s" % (
            bioawk_path, args.fna, 'fna.filter.bed')
        # cmd_info = "bash " + script_path + "filter_gtf_by_range.sh %s %s %s %s" % (
        #     bedtools_path, 'fna.filter.bed',
        #     args.input,
        #     'filter_range')
        os.system(cmd_bed)
        # os.system(cmd_info)
        cher2range=dict()
        with open('fna.filter.bed') as bed_r:
            for line in bed_r.readlines():
                line = line.strip().split('\t')
                cher2range[line[0]] = int(line[2])
        with open('filter_range', 'w') as filter_w, open(args.input, 'r') as in_r:
            scaf2line = defaultdict(dict)
            for line in in_r:
                if not line.startswith('#'):
                    tmp = line.strip().split('\t')
                    if len(tmp)>4 and tmp[0] in cher2range:
                        if int(tmp[4]) <= cher2range[tmp[0]]:
                            if not tmp[0] in scaf2line:
                                scaf2line[tmp[0]] = dict()
                            if not int(tmp[4]) in scaf2line[tmp[0]]:
                                scaf2line[tmp[0]][int(tmp[4])] = list()
                            if line not in scaf2line[tmp[0]][int(tmp[4])]:
                                scaf2line[tmp[0]][int(tmp[4])].append(line)
                            # filter_w.write(line)
            for scaf in sorted(scaf2line.keys()):
                for loc in sorted(scaf2line[scaf].keys()):
                    for i in scaf2line[scaf][loc]:
                        filter_w.write(i)
        if os.path.getsize('filter_range') == 0:
            raise Exception("参考基因组和gtf/gff文件不匹配，请核实后再投递！")
        else:
            ROCK = rock_index('filter_range', args.type, args.fna)
    else:
        ROCK = rock_index(args.input, args.type, args.fna)
    ROCK.run()



