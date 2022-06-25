# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/5/22 18:50

import re, os, Bio, argparse, sys, fileinput, urllib2
from Bio import SeqIO
import subprocess
import pickle
from biocluster.config import Config
from biocluster.agent import PickleConfig


def get_bed_from_gtf(bed_path, gtf_path):
    # todo 软件设置目录
    config = Config()
    os.environ['PATH'] = os.environ['PATH'] + '%s/bioinfo/align/bedops/bin' % config.SOFTWARE_DIR
    cmd = 'gtf2bed < %s > %s ' % (gtf_path, bed_path)
    try:
        subprocess.call(cmd, shell=True)
    except Exception as e:
        raise Exception('转换gtf文件到bed文件失败，命令为:\n %s' % cmd)


def get_fasta_index(fasta):
    config = Config()
    os.environ['PATH'] = os.environ['PATH'] + '%s/miniconda2/bin/' % config.SOFTWARE_DIR
    cmd = 'samtools faidx  %s ' % fasta
    try:
        subprocess.call(cmd, shell=True)
    except Exception as e:
        raise Exception('对fasta文件建索引失败，命令为:\n %s' % cmd)


def get_tbi(table_file, file_type):
    config = Config()
    os.environ['PATH'] = os.environ['PATH'] + '%s/bioinfo/align/htslib/bin/' % config.SOFTWARE_DIR
    cmd = 'tabix -p %s -f  %s ' % (file_type, table_file)
    try:
        subprocess.call(cmd, shell=True)
    except Exception as e:
        raise Exception('给文件建立tabix索引失败，命令为:\n %s' % cmd)


def get_bgzip_for_sorted_table_file(table_file):
    config = Config()
    os.environ['PATH'] = os.environ['PATH'] + '%s/bioinfo/align/htslib/bin/' % config.SOFTWARE_DIR
    cmd = 'bgzip -f  %s ' % (table_file)
    try:
        subprocess.call(cmd, shell=True)
    except Exception as e:
        raise Exception('bgzip 压缩文件失败，命令为:\n %s' % cmd)


def sort_vcf(vcf, sorted_vcf):
    tmp_head = vcf + ".head.temp"
    tmp_body = vcf + ".body.temp"
    tmp_sorted_body = vcf + ".sorted.body.temp"
    cmd = ''
    try:
        subprocess.call("grep '^#' %s > %s" % (vcf, tmp_head), shell=True)
        subprocess.call("grep -v  '^#' %s > %s" % (vcf, tmp_body), shell=True)
        cmd = 'sort -k 1,1 -k 2,2n  %s  > %s' % (tmp_body, tmp_sorted_body)
        subprocess.call(cmd, shell=True)
        subprocess.call('cat %s %s > %s' % (tmp_head, tmp_sorted_body, sorted_vcf), shell=True)
    except Exception as e:
        raise Exception('对vcf文件排序失败，命令为:\n %s' % cmd)


def sort_bed(bed, sorted_bed):
    tmp_head = bed + ".head.temp"
    tmp_body = bed + ".body.temp"
    tmp_sorted_body = bed + ".sorted.body.temp"
    cmd = ''
    try:
        subprocess.call("grep '^#' %s > %s" % (bed, tmp_head), shell=True)
        subprocess.call("grep -v  '^#' %s > %s" % (bed, tmp_body), shell=True)
        cmd = 'sort -k 1,1 -k 2,2n -k 3,3n  %s  > %s' % (tmp_body, tmp_sorted_body)
        subprocess.call(cmd, shell=True)
        subprocess.call('cat %s %s > %s' % (tmp_head, tmp_sorted_body, sorted_bed), shell=True)
    except Exception as e:
        raise Exception('对bed文排序失败，命令为:\n %s' % cmd)


def sort_gtf(gtf, sorted_gtf):
    tmp_head = gtf + ".head.temp"
    tmp_body = gtf + ".body.temp"
    tmp_sorted_body = gtf + ".sorted.body.temp"
    cmd = ''
    try:
        subprocess.call("grep '^#' %s > %s" % (gtf, tmp_head), shell=True)
        subprocess.call("grep -v  '^#' %s > %s" % (gtf, tmp_body), shell=True)
        cmd = 'sort -k 1,1 -k 4,4n -k 5,5n  %s  > %s' % (tmp_body, tmp_sorted_body)
        subprocess.call(cmd, shell=True)
        subprocess.call('cat %s %s > %s' % (tmp_head, tmp_sorted_body, sorted_gtf), shell=True)
    except Exception as e:
        raise Exception('对gtf文件排序失败，命令为:\n %s' % cmd)


def sort_bam(bam, sorted_bam):
    config = Config()
    os.environ['PATH'] = os.environ['PATH'] + '%s/miniconda2/bin/' % config.SOFTWARE_DIR
    cmd = 'samtools sort -o %s  %s ' % (sorted_bam, bam)
    try:
        subprocess.call(cmd, shell=True)
    except Exception as e:
        raise Exception('对bam文件排序失败，命令为:\n %s' % cmd)


def index_bam(bam):
    config = Config()
    os.environ['PATH'] = os.environ['PATH'] + '%s/miniconda2/bin/' % config.SOFTWARE_DIR
    cmd = 'samtools index  %s ' % bam
    try:
        subprocess.call(cmd, shell=True)
    except Exception as e:
        raise Exception('对bam文件建索引失败，命令为:\n %s' % cmd)


def get_fasta_chrom_size(fasta, chrom_size_file):
    temp_f = chrom_size_file + ".tmp"
    fw = open(temp_f, 'w')
    for seq in SeqIO.parse(fasta, 'fasta'):
        fw.write('%s\t%s\n' % (seq.id, len(seq.seq)))
    fw.close()
    subprocess.call('sort -k 1,1 %s > %s' % (temp_f, chrom_size_file), shell=True)
    subprocess.call('rm -rf %s' % temp_f)


def get_bigbed_from_bed(bed, chrom_size, bigbed):
    '''
    
    需要先对bed文件进行排序, 还需要chrom.size文件
    :param bed:
    :param bigbed:
    :return:
    '''
    config = Config()
    os.environ['PATH'] = os.environ['PATH'] + '%s/bioinfo/align/ucsc_tools/' % config.SOFTWARE_DIR
    cmd = 'bedToBigBed  %s  %s  %s' % (bed, chrom_size, bigbed)
    try:
        subprocess.call(cmd, shell=True)
    except Exception as e:
        raise Exception('bed转bigbed文件失败，命令为:\n %s' % cmd)


def get_wig_from_bam(bam, out_dir):
    '''
    1. 需要设置好samtools路径放在环境变量里
    2. 需要设置好ucsc的wigToBigWig工具所在文件夹在PATH路径里
    3. 输入的bam必须是sort后的，而且index过的， 而且bam的bai文件要放在与自己同一个文件夹里
    4. 使用的是https://github.com/MikeAxtell/bam2wig/blob/master/bam2wig bam2wig的脚本
    
    :param bam:
    :param out_dir:
    :return:
    '''
    
    config = Config()
    os.environ['PATH'] = os.environ['PATH'] + '%s/bioinfo/align/ucsc_tools/' % config.SOFTWARE_DIR
    cmd = 'wigToBigWig -D  %s  %s' % (out_dir, bam)
    try:
        subprocess.call(cmd, shell=True)
    except Exception as e:
        raise Exception('bam转wig文件失败，命令为:\n %s' % cmd)


def get_fasta_chrom_sizes(fasta, chrom_sizes):
    temp_f = chrom_sizes + ".tmp"
    fw = open(temp_f, 'w')
    for seq in SeqIO.parse(fasta, 'fasta'):
        fw.write('%s\t%s\n' % (seq.id, len(seq.seq)))
    fw.close()
    subprocess.call('sort -k 1,1 %s > %s' % (temp_f, chrom_sizes), shell=True)
    subprocess.call('rm -rf %s' % temp_f)


def get_init_locus(ref_fasta):
    seq_records = dict([(seq_record.id, len(seq_record.seq)) for seq_record in SeqIO.parse(ref_fasta, "fasta")])
    max_len = max(seq_records.values())
    max_id = ''
    for seq_id in seq_records:
        if seq_records[seq_id] == max_len:
            max_id = seq_id
            break
    start = round(int(max_len) / 3)
    end = round(3 * int(max_len) / 7)
    return '%s:%s-%s' % (max_id, start, end)


