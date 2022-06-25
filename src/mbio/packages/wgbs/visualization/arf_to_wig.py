# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import re
import os
import sys
import collections
import numpy as np
import regex
import csv
from Bio import SeqIO
from reportlab.graphics import renderPDF, renderPM
from svglib.svglib import svg2rlg


class ArfToWig(object):
    def __init__(self, sample, types=None):
        self.arf = ""
        self.wig = ""
        self.step_width = 5000
        self.mapping = collections.OrderedDict()
        self.fa = ""
        self.sample = sample
        self.types = types
        self.indexs = []
        self.chr_list = []

    def set_arf(self, arf):
        '''
        设置arf路径
        '''
        self.arf = arf

    def set_fa(self, fa):
        '''
        设置fa路径
        '''
        self.fa = fa

    def init_mapping(self):
        '''
        初始化列表
        '''
        index = 0
        for seq in SeqIO.parse(self.fa, "fasta"):
            chrom = seq.id
            length = len(seq.seq)
            start = 1
            while True:
                end = start + self.step_width - 1
                if end > length:
                    end = length
                    pos = (start + end)/2
                    self.mapping.update({chrom + "_" + str(start): {
                        "chr": chrom,
                        "start": start,
                        "end": end,
                        "pos": pos,
                        "num_pos": [],
                        "num_neg": []
                    }})
                    break
                else:
                    pos = (start + end)/2
                    self.mapping.update({chrom + "_" + str(start): {
                        "chr": chrom,
                        "start": start,
                        "end": end,
                        "pos": pos,
                        "num_pos": [],
                        "num_neg": []
                    }})
                    start = end + 1

    def parser_arf(self):
        '''
        统计arf文件
        '''
        with open(self.arf, 'r') as f:
            for dic in csv.DictReader(f, delimiter="\t"):
                # print dic
                if self.types != None:
                    if dic['c_content'] != self.types:
                        continue

                if type(sample) == tuple:
                    rate = np.mean([float(dic[s]) for s in sample[1]])
                else:
                    rate = float(dic[self.sample])
                start = int(dic['position'])
                end = int(dic['position'])
                chrom = dic['chromosome']

                start_pre = start/self.step_width * self.step_width + 1
                index_pre = chrom + '_' + str(start_pre)
                index_next = chrom + '_' + str(start_pre + self.step_width)
                if dic["strand"] == "+":
                    num_strand = "num_pos"
                else:
                    num_strand = "num_neg"

                if num_strand in self.mapping[index_pre]:
                    self.mapping[index_pre][num_strand].append(rate)
                else:
                    self.mapping[index_pre][num_strand] = [rate]

    def write_file(self, stat_file):
        '''
        输出统计相关结果
        '''
        with open(stat_file, 'wb') as f:
            for index,values in self.mapping.items():
                f.write("\t".join([values['chr'], str(values['start']), str(values['end']), str(values['pos']), str(values['num_pos']), str(values['num_neg'])]) + '\n')

    def write_window(self):
        '''
        输出统计相关结果
        '''
        if type(self.sample) == tuple:
            sample = self.sample[0]
        else:
            sample = self.sample
        with open(sample + '.pos.window', 'wb') as f_pos, open(sample + '.neg.window', 'wb') as f_neg:
            pos_list = [value_dict['num_pos'] for value_dict in self.mapping.values()]
            neg_list = [value_dict['num_neg'] for value_dict in self.mapping.values()]
            num_max = max(pos_list + neg_list)
            chr_set = set(self.chr_list)
            for index,values in self.mapping.items():
                if values['chr'] in chr_set:
                    if len(values['num_pos']) > 0:
                        pos = np.mean(values['num_pos'])
                    else:
                        pos = 0
                    if len(values['num_neg']) > 0:
                        neg = np.mean(values['num_neg'])
                    else:
                        neg = 0
                    f_pos.write("\t".join([values['chr'], str(values['start']), str(values['end'] - 1), str(pos), "fill_color=col1"])  + '\n')
                    f_neg.write("\t".join([values['chr'], str(values['start']), str(values['end'] - 1), str(neg), "fill_color=col3"])  + '\n')

class gtf2label(object):
    def __init__(self, gtf):
        self.gtf = gtf

    def convert2label(self, level="gene", show="name"):
        with open(self.gtf, "r") as f, open("label.txt", "w") as fo:
            fo.write("#chr\tstart\tend\tlabel\n")
            for line in f:
                cols = line.strip("\n").split("\t")
                if cols[2] == "gene":
                    txpt_id = ''
                    gene_id = ''
                    gene_name = ""
                    print cols[8]
                    content_m = regex.match(
                        r'^(.*;)*\s*((gene_id|gene_name)\s+?\"(\S+?)\");.*((gene_name|gene_id)\s+?\"(\S+?)\");(.*;)*$',
                        cols[8])
                    if content_m:
                        if 'gene_name' in content_m.captures(6):
                            gene_name = content_m.captures(7)[0]
                        elif 'gene_name' in content_m.captures(3):
                            gene_name = content_m.captures(4)[0]
                        else:
                            gene_name = ""
                        fo.write("{}\t{}\t{}\t{}\n".format(cols[0], cols[3], cols[4], gene_name))
                    else:
                        pass

class bed2label(object):
    def __init__(self, bed, step_width):
        self.bed = bed
        self.step_wigth = step_width
        self.genome = collections.OrderedDict()
        self.region = dict()
        self.get_bed_info()

    def get_bed_info(self):
        with open(self.bed, 'r') as f:
            for line in f:
                cols = line.strip().split("\t")
                if cols[0] in self.region:
                    self.region[cols[0]].append(cols)
                else:
                    self.region[cols[0]] = [cols]

    def convert2band(self, fai, chr_band):
        with open(fai, 'r') as f, open(chr_band, 'w') as fo:
            chr_num = 0
            for line in f:
                chr_num += 1
                band_num = 0
                cols = line.strip().split("\t")
                self.genome[cols[0]] = int(cols[1])
                if cols[0] in self.region:
                    fo.write("\t".join(["chr", "-", cols[0], cols[0], "0", cols[1], cols[0]]) + "\n")
                    for line in self.region[cols[0]]:
                        band_num += 1
                        start = int(line[1])/int(self.step_wigth) * int(self.step_wigth)  + 1
                        end = int(line[2])/int(self.step_wigth) * int(self.step_wigth) + int(self.step_wigth) - 1
                        band_col = str(chr_num % 6 + 1)
                        fo.write("\t".join(["band",
                                            cols[0],
                                            "band{}".format(band_num),
                                            "band{}".format(band_num),
                                            str(start),
                                            str(end),
                                            "col{}".format(band_col)
                        ])  + "\n")


    def convert2label(self, level="gene", show="name"):
        with open(self.bed, "r") as f, open("label.txt", "w") as fo:
            fo.write("#chr\tstart\tend\tlabel\n")
            for line in f:
                cols = line.strip("\n").split("\t")
                fo.write("{}\t{}\t{}\t{}\n".format(cols[0], cols[1], cols[2], cols[3]))

def svg_convert(fin, fo):
    '''
    修改图片转换为并行
    '''

    drawing = svg2rlg(fin)
    renderPDF.drawToFile(drawing, fo)

if __name__ == "__main__":
    test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/wgbs/'
    '''
    G = bed2label(test_dir + 'rCRS.GRCh38.96.bed', 10)
    G.convert2band(test_dir + 'rCRS.fasta.fai', "chr.band.txt")
    G.convert2label()
    '''
    for sample in [("C", ["C1", "C2", "C3", "C4", "C5"]), ("L", ["L1", "L2", "L3", "L4", "L5", "L6"])]:
        # sample = "C1"
        A = ArfToWig(sample)
        A.step_width = 10
        A.set_fa(test_dir + 'rCRS.fasta')
        A.set_arf(test_dir + 'mth_ratio_gene.txt')
        A.init_mapping()
        # print A.mapping
        A.parser_arf()
        A.chr_list=['rCRS']
        A.write_window()

        circos_cmd = "/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/perl /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/wgbs/visualization/draw_circos.pl --windows 10 --pos {}.pos.window --neg {}.neg.window --label label.txt --chrlist rCRS.fasta.fai --chr_band chr.band.txt --outdir {}".format(sample[0], sample[0], sample[0])
        print circos_cmd
        os.system(circos_cmd)
        convert_cmd = "rsvg-convert -f pdf -o {}.circos.pdf {}/circos.svg".format(sample[0], sample[0])
        os.system(convert_cmd)
        # rsvg-convert -f pdf -o C1/circos.pdf C1/circos.svg
        # svg_convert("{}/circos.svg".format(sample), "{}.pdf".format(sample))

        # sample = "C1"
        A = ArfToWig(sample, 'CG')
        A.step_width = 10
        A.set_fa(test_dir + 'rCRS.fasta')
        A.set_arf(test_dir + 'mth_ratio_gene.txt')
        A.init_mapping()
        # print A.mapping
        A.parser_arf()
        A.chr_list=['rCRS']
        A.write_window()


        circos_cmd = "/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/perl /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/wgbs/visualization/draw_circos.pl --windows 10 --pos {}.pos.window --neg {}.neg.window --label label.txt --chrlist rCRS.fasta.fai --chr_band chr.band.txt --outdir {}".format(sample[0], sample[0], sample[0] + "_cg")
        print circos_cmd
        os.system(circos_cmd)

        convert_cmd = "rsvg-convert -f pdf -o {}.circos.pdf {}/circos.svg".format(sample[0] + "_cg", sample[0] + "_cg")
        os.system(convert_cmd)
        # svg_convert("{}/circos.svg".format(sample + "_cg"), "{}.pdf".format(sample + "_cg"))

        # G = gtf2label(test_dir + 'rCRS.GRCh38.96.gtf')
        # G.convert2label()
    '''
    for sample in ["C1", "C2", "C3", "C4", "C5", "L1", "L2", "L3", "L4", "L5", "L6"]:
        # sample = "C1"
        A = ArfToWig(sample)
        A.step_width = 10
        A.set_fa(test_dir + 'rCRS.fasta')
        A.set_arf(test_dir + 'mth_ratio_gene.txt')
        A.init_mapping()
        # print A.mapping
        A.parser_arf()
        A.chr_list=['rCRS']
        A.write_window()

        circos_cmd = "/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/perl /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/wgbs/visualization/draw_circos.pl --windows 10 --pos {}.pos.window --neg {}.neg.window --label label.txt --chrlist rCRS.fasta.fai --chr_band chr.band.txt --outdir {}".format(sample, sample, sample)
        print circos_cmd
        os.system(circos_cmd)
        convert_cmd = "rsvg-convert -f pdf -o {}.circos.pdf {}/circos.svg".format(sample, sample)
        os.system(convert_cmd)
        # rsvg-convert -f pdf -o C1/circos.pdf C1/circos.svg
        # svg_convert("{}/circos.svg".format(sample), "{}.pdf".format(sample))

        # sample = "C1"
        A = ArfToWig(sample, 'CG')
        A.step_width = 10
        A.set_fa(test_dir + 'rCRS.fasta')
        A.set_arf(test_dir + 'mth_ratio_gene.txt')
        A.init_mapping()
        # print A.mapping
        A.parser_arf()
        A.chr_list=['rCRS']
        A.write_window()


        circos_cmd = "/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/perl /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/wgbs/visualization/draw_circos.pl --windows 10 --pos {}.pos.window --neg {}.neg.window --label label.txt --chrlist rCRS.fasta.fai --chr_band chr.band.txt --outdir {}".format(sample, sample, sample + "_cg")
        print circos_cmd
        os.system(circos_cmd)

        convert_cmd = "rsvg-convert -f pdf -o {}.circos.pdf {}/circos.svg".format(sample + "_cg", sample + "_cg")
        os.system(convert_cmd)
        # svg_convert("{}/circos.svg".format(sample + "_cg"), "{}.pdf".format(sample + "_cg"))

        # G = gtf2label(test_dir + 'rCRS.GRCh38.96.gtf')
        # G.convert2label()
    '''
