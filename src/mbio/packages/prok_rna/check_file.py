# -*- coding: utf-8 -*-
# __author__ : 'shicaiping'
# __date__: 202108

from Bio import SeqIO
from BCBio import GFF
import os
import subprocess
from biocluster.config import Config
import logging


class CheckFile(object):
    """
    Used to check whether the file format of genome fasta and gtf/gff file meet the standard or not.
    """

    def __init__(self, genome, gff=None, gtf=None, work_dir=None):
        logging.basicConfig()
        self.logger = logging.getLogger()
        self.logger.setLevel(logging.INFO)
        self.config = Config()
        self.gff_read = self.config.SOFTWARE_DIR + '/bioinfo/rna/cufflinks-2.2.1/gffread'
        self.gtf2bed_path = self.config.SOFTWARE_DIR + "/bioinfo/rna/scripts/gtf2bed.pl"
        self.genome = genome
        if gff:
            self.anno_file = gff
            self.anno_type = "gff"
        else:
            self.anno_file = gtf
            self.anno_type = "gtf"
        if work_dir:
            self.work_dir = work_dir
        else:
            self.work_dir = os.getcwd()

    def run(self):
        info4 = self.check_genome()
        if not isinstance(info4, bool):
            return info4
        if self.anno_type == "gff":
            info = self.check_gff()
            if not isinstance(info, bool):
                return info
            info1 = self.to_gtf()
            if not isinstance(info1, bool):
                return info1
        info2 = self.check_gtf()
        if not isinstance(info2, bool):
            return info2
        info3 = self.to_bed()
        if not isinstance(info3, bool):
            return info3
        self.logger.info("完成文件检查！")
        return True

    def check_genome(self):
        self.logger.info("开始检查参考基因组文件和注释文件是否配套")
        chr_range = dict()
        genome_range = os.path.join(self.work_dir, "genome_range.txt")
        with open(genome_range, "w") as w:
            for seq_record in SeqIO.parse(self.genome, "fasta"):
                if seq_record.id not in chr_range:
                    chr_range[seq_record.id] = len(seq_record.seq)
                    w.write(seq_record.id + ":0.." + str(len(seq_record.seq)) + "\n")
        with open(self.anno_file, "r") as f:
            not_include = list()
            out_of_range = list()
            for line in f:
                if line.startswith("#") or not len(line.strip()):
                    continue
                if line.strip().split("\t")[0] not in chr_range:
                    if line.strip().split("\t")[0] not in not_include:
                        not_include.append(line.strip().split("\t")[0])
                elif int(line.strip().split("\t")[4]) > chr_range[line.strip().split("\t")[0]]:
                    if line.strip().split("\t")[0] not in out_of_range:
                        out_of_range.append(line.strip().split("\t")[0])
            if not_include and out_of_range:
                return '{}文件和{}不匹配, {}在参考基因组文件中不存在; 且{}超过参考基因组坐标范围！'.format(os.path.basename(self.genome), os.path.basename(self.anno_file),
                                                                         ",".join(not_include),
                                                                         ",".join(out_of_range))
            elif not_include:
                return '{}文件和{}不匹配, {}在参考基因组文件中不存在！'.format(os.path.basename(self.genome), os.path.basename(self.anno_file),
                                                                     ",".join(not_include))
            elif out_of_range:
                return '{}超过参考基因组坐标范围！'.format(",".join(out_of_range))
        return True

    def check_gff(self):
        self.logger.info("开始检查GFF文件格式是否符合规范")
        in_handle = open(self.anno_file)
        rec_list = []
        try:
            for rec in GFF.parse(in_handle):
                try:
                    rec_list.append(rec)
                except:
                    return "{}文件格式错误！".format(self.anno_file)
        except:
            return "{}文件格式错误！".format(self.anno_file)
        return True

    def to_gtf(self):
        self.logger.info("开始转换GFF文件为GTF文件")
        to_gtf_cmd = '%s %s -T -o %s  ' % (self.gff_read, self.anno_file, os.path.join(self.work_dir, 'ref.gtf'))
        try:
            subprocess.check_output(to_gtf_cmd, shell=True)
        except subprocess.CalledProcessError as exc:
            if exc.returncode is not None:
                if not isinstance(info, string):
                    return "GFF格式转换GTF失败，请检查格式是否符合规范！"
                else:
                    return exc.returncode
            else:
                return "GFF格式转换GTF失败，请检查格式是否符合规范！"
        return True

    def check_gtf(self):
        self.logger.info("开始检查GTF文件格式是否符合规范")
        flag1 = True
        flag2 = True
        flag3 = True
        if self.anno_type == "gtf":
            gtf = self.anno_file
        else:
            gtf = os.path.join(self.work_dir, 'ref.gtf')
        with open(gtf, 'r') as f:
            for line in f:
                if len(line.strip().split("\t")) == 9:
                    flag1 = False
                    if 'gene_id ' in line.strip().split("\t")[8]:
                        flag2 = False
                if len(line.strip().split("\t")) >= 3:
                    if line.strip().split("\t")[2].lower() == "cds":
                        flag3 = False
        if flag1:
            return "{}文件格式不对，至少包含9列信息！".format(self.anno_file)
        if flag2:
            return "{}文件格式不对，attribute列至少包含gene_id信息".format(self.anno_file)
        if flag3:
            return "{}文件格式不对，未包含CSD信息".format(self.anno_file)
        return True

    def to_bed(self):
        self.logger.info("开始转换GTF文件为BED文件")
        if self.anno_type == "gtf":
            gtf_path = self.anno_file
        else:
            gtf_path = os.path.join(self.work_dir, 'ref.gtf')
        bed_path = os.path.join(self.work_dir, 'ref.bed')
        cmd = "perl {} {} > {}".format(self.gtf2bed_path, gtf_path, bed_path)
        try:
            subprocess.check_output(cmd, shell=True)
        except:
            return "GFF格式转换为BED失败，请检查格式是否符合规范！"
        return True


if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description='This script is used to check whether the file format of genome fasta and gtf/gff file meet the standard or not.')
    parser.add_argument('-genome', type=str, required=True, help='full path of genome fasta file.')
    parser.add_argument('-gff', type=str, help='full path of gtf file.')
    parser.add_argument('-gtf', type=str, help='full path of gff file. When gff is not provided.')
    parser.add_argument('-work_dir', type=str, help='full path of directory where temporary files are stored.')
    args = parser.parse_args()
    if args.gff:
        check = CheckFile(args.genome, gff=args.gff, work_dir=args.work_dir)
    elif args.gtf:
        check = CheckFile(args.genome, gtf=args.gtf, work_dir=args.work_dir)
    else:
        raise Exception("Must be provide at least one of gtf/gff.")
    check.run()
