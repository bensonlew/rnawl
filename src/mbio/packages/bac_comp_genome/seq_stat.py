# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os
import argparse
from Bio import SeqIO

"""
说明：
本脚本处理complete和draft两种情况的统计；
如果为complete情况，分别统计每条序列的信息；
如果为draft的情况，需要将contig和scaffold分别进行统计，最后将结果统计到一起，整理为一行
一个基因组为一行数据，最后统计出所有基因组的数据
"""


def seq_stat(input, seq, raw, output):
    """
    输入参数文件，seq序列的文件夹，序列原始序列的文件夹，输出文件夹
    """
    outw = open(output, 'w')
    outw.write("Assembly_accession\tTotal_length\tGenome_Name\tGenome_Type\tGenome_length\tGC\tChrom_num\tPlasmid_num\tScaffold_num\tContig_num\tScaffold_n50\tContig_n50\tSeq_status\tSpecies\tGene_prefix\n")
    with open(input, 'r') as f:
        lines = f.readlines()
        for line in lines[1:]:
            line = line.strip().split("\t")
            type = line[1]
            assembly_type = covert_dict(line[2])
            gene_prefix = line[-1]
            species = line[-2]
            if type in ['seq', 'gff']:
                if assembly_type in ['draft']:
                    assembly_name = line[0]
                    genome_name = line[3]
                    genome_type = line[4]
                    seq_path = os.path.join(raw, genome_name + ".fna")
                    try:
                        (seq_num, gc_content, genome_length) = draft_stat(seq_path)
                        (scaffold_num,scaffold_50) = cal_scaffold(genome_length,seq_path)
                        (total_contig_num, contig_n50) = cal_contig(genome_length,seq_path)
                    except:
                        seq_path = line[5]
                        (seq_num, gc_content, genome_length) = draft_stat(seq_path)
                        (scaffold_num,scaffold_50) = cal_scaffold(genome_length,seq_path)
                        (total_contig_num, contig_n50) = cal_contig(genome_length,seq_path)
                    outw.write(assembly_name +"\t"+ str(genome_length) + "\t"+ str(genome_name)+ "\t"+ str(genome_type)+ "\t"+ str(genome_length) +"\t"+ str(gc_content) +"\t"+ "-\t-\t" + str(scaffold_num) +"\t"+ str(total_contig_num) +"\t"+ str(scaffold_50) +"\t"+ str(contig_n50) +"\tdraft\t"+ species +"\t"+ gene_prefix +"\n")
                elif assembly_type in ['complete', 'chromosome']:
                    assembly_name = line[0]
                    genome_str = line[3]
                    genome_type = line[4]
                    genome_type_list = genome_type.split(",")
                    genome_list = genome_str.split(",")
                    genome_type_list_dict = zip(genome_list, genome_type_list)
                    total_length = 0
                    chromo_num = 0
                    plas_num = 0
                    genome_accession = []
                    genome_length_list = []
                    gc_list = []
                    for genome in genome_list:
                        list_index = genome_list.index(genome)
                        new_genome_type = genome_type_list_dict[list_index][1]
                        if new_genome_type in ['chromosome', 'Chromosome', 'CHROMOSOME', 'chrom', 'Chrom']:
                            chromo_num += 1
                        elif new_genome_type in ['Plasmid', 'plasmid', 'PLASMID']:
                            plas_num += 1
                        genome_name = genome
                        assembly_path = os.path.join(seq, assembly_name)
                        try:
                            seq_path = os.path.join(assembly_path, genome + ".fna")
                            seq_num, gc_content, genome_length = draft_stat(seq_path)
                        except:
                            index = genome_list.index(genome)
                            seq_path_list = line[5].split(",")
                            seq_path = seq_path_list[index]
                            print(seq_path)
                            seq_num, gc_content, genome_length = draft_stat(seq_path)
                        genome_accession.append(genome_name)
                        genome_length_list.append(str(genome_length))
                        gc_list.append(str(gc_content))
                        total_length += genome_length
                    print("genome_length_list: %s" %genome_length_list)
                    print("genome_list: %s" %genome_list)
                    print("genome_type_list: %s" %genome_type_list)
                    print("gc_list: %s" %gc_list)
                    outw.write(assembly_name +"\t"+ str(total_length) +"\t"+ ','.join(genome_list) +"\t"+ ','.join(genome_type_list) +"\t"+ ','.join(genome_length_list) +"\t"+ ','.join(gc_list) +"\t"+ str(chromo_num) +"\t"+ str(plas_num) +"\t-\t-\t-\t-\t"+ assembly_type +"\t"+ species +"\t"+ gene_prefix +"\n")
                else:
                    raise Exception('基因组的类型设置不正确')
            else:
                pass


def cal_scaffold(genome_length, seq_path):
    """
    计算scaffold的N50
    :param genome_length: 基因组大小
    :param seq_path: 基因组序列路径
    :return:
    """
    scaffold_length = genome_length * 0.5
    seq_list = []
    scaf_num = 0
    scaffold_n50 = 0
    scaffold_num = 0
    for seq_record in SeqIO.parse(seq_path, "fasta"):
        seq_list.append(seq_record.seq)
        scaffold_num += 1
    seq_list.sort(key=lambda i:len(i), reverse=True)
    for new_seq in seq_list:
        seq_length = len(new_seq)
        scaf_num += seq_length
        if scaf_num < scaffold_length:
            continue
        else:
            scaffold_n50 = seq_length
            break
    return(scaffold_num,scaffold_n50)


def cal_contig(genome_length, seq_path):
    """
    计算contig的N50
    :param genome_length: 基因组大小
    :param seq_path: 基因组序列路径
    :return:
    """
    contig_length = genome_length * 0.5
    seq_contig_list = []
    contig_num = 0
    contig_n50 = 0
    total_contig_num = 0
    for seq_record in SeqIO.parse(seq_path, "fasta"):
        seq = str(seq_record.seq)
        replace_seq = seq.replace("n", "N")
        seq_list = replace_seq.split("N")
        for ss in seq_list:
            # if ss not in seq_contig_list:
            seq_contig_list.append(ss)
            total_contig_num += 1
    seq_contig_list.sort(key=lambda i:len(i), reverse=True)
    for new_seq in seq_contig_list:
        seq_length = len(new_seq)
        contig_num += seq_length
        if contig_num < contig_length:
            continue
        else:
            contig_n50 = seq_length
            break
    return(total_contig_num,contig_n50)


def draft_stat(seq_path):
    """
    统计序列文件的大小、gc、序列数
    :param seq_path: fasta序列的路径
    :return:
    """
    genome_length = 0
    seq_num = 0
    gc_num = 0
    for seq_record in SeqIO.parse(seq_path, "fasta"):
        seq = seq_record.seq
        G_num = seq.count("G")
        g_num = seq.count("g")
        C_num = seq.count("C")
        c_num = seq.count("c")
        gc_num = gc_num + G_num + C_num + g_num + c_num
        genome_length += len(seq)
        seq_num += 1
    gc_content = float(gc_num) / genome_length
    return(seq_num, gc_content, genome_length)


def covert_dict(key):
    """
    对输入的key进行转换
    :param key: 输入类型
    :return:
    """
    if key in ["Draft", "draft", "DRAFT"]:
        assembly_type = "draft"
    elif key in ["Complete Genome", "complete", "Complete", "COMPLETE"]:
        assembly_type = "complete"
    elif key in ["Chromosome", "chromosome"]:
        assembly_type = "chromosome"
    else:
        assembly_type = ""
    return assembly_type


def main():
    parse = argparse.ArgumentParser()
    parse.add_argument('-i', metavar='[total_genome]',required=True, help='input params file')
    parse.add_argument('-seq', metavar='[seq dir]', required=True, help='input seq dir')
    parse.add_argument('-total', metavar='[total_raw_dir]', required=True, help='input total raw dir')
    parse.add_argument('-o', metavar='[outfile]', required=True, help='input outfile')
    args = parse.parse_args()
    input = args.i
    seq = args.seq
    raw = args.total
    output = args.o
    seq_stat(input, seq, raw, output)

if __name__ == '__main__':
    main()
