# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue,shicaiping'

import fileinput
import os
import re
import time
from collections import defaultdict

import regex
from Bio import SeqIO


def step_count(fasta_file, fasta_to_txt, group_num, step, stat_out):
    """
    步长统计
    :param fasta_file: 输入的fa文件
    :param fasta_to_txt:输出的统计数据txt
    :param group_num:按照步长统计几组
    :param step:步长
    :param stat_out:统计的数据汇总信息txt（fasta_to_txt文件的汇总）
    :return:
    """
    if not (os.path.isfile(fasta_to_txt) and time.time() - os.stat(fasta_to_txt).st_mtime < 3600):
        with open(fasta_to_txt, "w") as f:
            for seq_record in SeqIO.parse(fasta_file, "fasta"):
                f.write("{}\t{}\n".format(seq_record.description.strip().split(" ")[0], len(seq_record)))
    with open(fasta_to_txt, "r") as r, open(stat_out, "a") as w:
        sample_name = os.path.basename(fasta_file).split('.fa')[0]
        w.write("{}\n".format(sample_name))
        trans_list = []
        amount_group = []
        element_set = set("")
        for line in r:
            line = line.strip().split("\t")
            number = line[1]
            trans_list.append(number)
        for f in trans_list:
            for i in range(group_num):
                if (int(f) >= (i * step)) and (int(f) < ((i + 1) * step)):
                    amount_group.append(i)
                element_set.add(i)
        amount_group.sort()
        top_sum = 0
        for i in element_set:
            num_statistics = amount_group.count(i)
            if str(i) == '0':
                area_line = str(i * step) + "~" + str((i + 1) * step) + "\t" + str(num_statistics) + "\n"
                w.write(area_line)
                top_sum += int(num_statistics)
            elif i < (group_num - 1):
                area_line = str(i * step + 1) + "~" + str((i + 1) * step) + "\t" + str(num_statistics) + "\n"
                w.write(area_line)
                top_sum += int(num_statistics)
            else:
                area_line = ">" + str(i * step) + "\t" + str(len(trans_list) - int(top_sum)) + "\n"
                end_line = "total" + "\t" + str(len(trans_list)) + "\n"
                w.write(area_line)
                w.write(end_line)


def merged_add_code(trans_file, tmap_file, new_trans):
    def get_info_dic_from_tmap(srcfile):  # 转录本信息
        dic = dict()
        for line in fileinput.input(srcfile):
            arr = line.strip().split("\t")
            key = arr[4]  # transctipt id
            value = arr[2]
            dic[key] = value
        return dic

    candidateDic = get_info_dic_from_tmap(tmap_file)
    p = re.compile(r'transcript_id')
    fw = open(new_trans, "w+")
    for line in fileinput.input(trans_file):
        m = re.match("#.*", line)
        if not m:
            line1 = line.strip()
            arr = line1.split("\t")
            description = arr[8]
            desc_array = description.strip().split(";")
            for tag in desc_array:
                tagpair = tag.strip().split("\"")
                if len(tagpair) == 3:
                    if re.match(p, tagpair[0]):
                        if candidateDic.has_key(tagpair[1]):
                            newline1 = line1 + " class_code \"" + candidateDic[tagpair[1]] + "\";\n"
                            fw.write(newline1)

    fw.close()


def count_trans_or_exons(input_file, step, count_file, final_files):
    """
    统计对应关系的信息
    :param input_file:第一列：基因id，第二列：长度
    :param step:步长，范围内的求和
    :param count_file:第一列：基因（转录本）对应的转录本（外显子）的个数，第二列：这样的基因（转录本）有多少个，第三列：list,分别是哪些基因（转录本）id
    :param final_files:合并之后的count_file
    """
    with open(input_file, "r") as r:
        lines = r.readlines()
        if len(lines) >= 3:
            f1 = open(count_file, 'w')
            f2 = open(final_files, 'w')
            dic = {}
            group = 0
            list_set = set()
            for line in fileinput.input(input_file):
                lines = line.strip().split("\t")
                dic[lines[0]] = lines[1]
                list_set.add(int(lines[1]))
            sort_set = list(list_set)
            for i in sorted(sort_set):
                value = 0
                ids_list = []
                for key in dic.keys():
                    if int(dic[key]) == i:
                        value += 1
                        ids_list.append(key)
                    else:
                        pass
                if int(step) == 1:
                    new_line = str(i) + "\t" + str(value) + "\t" + str(ids_list) + "\n"
                    f1.write(new_line)
                    f2.write(new_line)
                else:
                    new_line = str(i) + "\t" + str(value) + "\t" + str(ids_list) + "\n"
                    f1.write(new_line)
            f1.close()
            if int(step) != 1:
                max_num = sorted(sort_set)[-1]
                group = max_num / step
                mod = max_num % step
                if mod != 0:
                    group += 1
                else:
                    pass
            for n in range(group):
                f3 = open(count_file, 'r+')
                final_value = 0
                final_ids = []
                for line in f3:
                    num = line.strip().split("\t")[0]
                    value = line.strip().split("\t")[1]
                    id_s = line.strip().split("\t")[2]
                    ids_list = id_s.strip("[").strip("]").strip("'").split(",")
                    if (int(num) > (n * step)) and (int(num) <= ((n + 1) * step)):
                        final_value += int(value)
                        for id_num in ids_list:
                            final_ids.append(id_num)
                if n == 0:
                    area_line = str(n * step) + "~" + str((n + 1) * step) + "\t" + str(final_value) + "\t" + str(
                        final_ids) + "\n"
                else:
                    area_line = str(n * step + 1) + "~" + str((n + 1) * step) + "\t" + str(final_value) + "\t" + str(
                        final_ids) + "\n"
                f2.write(area_line)
                f3.close()
            f2.close()


def gene_trans_exon(merged_gtf, method, gene_trans_file, trans_exon_file):
    fw1 = open(gene_trans_file, 'w+')
    fw2 = open(trans_exon_file, 'w+')
    gene_set_dic = defaultdict(set)
    gene_dic = defaultdict(int)
    trans_dic = defaultdict(int)
    for line in fileinput.input(merged_gtf):
        m = regex.search(r'^[^#]\S*\t(\S+\t){7}.*?gene_id\s+\"(\S+)\";.*\s+transcript_id\s+\"(\S+)\";*', line)
        if m:
            seq_type = m.captures(1)[1]
            txpt_id = m.captures(3)[0]
            gene_id = m.captures(2)[0]
            if str(method).lower() == "stringtie":
                if seq_type.split("\t")[0] == "exon":
                    trans_dic[txpt_id] += 1
            elif str(method).lower() == "cufflinks":
                trans_dic[txpt_id] += 1
            gene_set_dic[gene_id].add(txpt_id)
    for key in gene_set_dic.keys():
        gene_dic[key] = len(gene_set_dic[key])
        gene_line = key + "\t" + str(gene_dic[key]) + "\n"
        fw1.write(gene_line)
    for key in trans_dic.keys():
        trans_line = key + "\t" + str(trans_dic[key]) + "\n"
        fw2.write(trans_line)
    fw1.close()
    fw2.close()


def class_code_count(gtf_file, code_num_trans):
    f = open(gtf_file, "rb")
    txpt_cls_content = []
    cls_content = set()
    for line in f:
        m = re.match("#.*", line)
        if not m:
            nine_line = line.strip().split("\t")[-1]
            n = re.search(r'\s*transcript_id\s+\"(\S+)\";.*\s*class_code\s+\"(\S+)\";*', nine_line)
            if n:
                cls_content.add(n.group(2))
                new_line = 'transcript_id "' + n.group(1) + '";\t class_code "' + n.group(2) + '"\n'
                txpt_cls_content.append(new_line)

    cls_txpt_set_dic = {}
    for cls in cls_content:
        cls = cls.strip()
        if cls:
            cls_txpt_set_dic[cls] = {'txpt_set': set(), 'count': 0}

    for record in txpt_cls_content:
        m = re.search(r'\s*transcript_id\s+\"(\S+)\";\s*class_code\s+\"(\S+)\"', record.strip())
        if m:
            cls = m.group(2)
            txpt = m.group(1)
            cls_txpt_set_dic[cls]['txpt_set'].add(txpt)

    fw = open(code_num_trans, 'wb')
    for cls in cls_txpt_set_dic.keys():
        cls_txpt_set_dic[cls]['count'] = len(cls_txpt_set_dic[cls]['txpt_set'])
        newline = '{}\t{}\t{}\n'.format(cls, ','.join(cls_txpt_set_dic[cls]['txpt_set']),
                                        str(cls_txpt_set_dic[cls]['count']))
        fw.write(newline)
    fw.close()


def tidy_code_num(all_transcripts_fa, raw_code_num_txt, new_code_num_txt):
    transcript_id_set = {r.id for r in SeqIO.parse(all_transcripts_fa, 'fasta')}
    transcript_id_to_code = dict()
    for line in open(raw_code_num_txt):
        eles = line.strip().split('\t')
        if len(eles) == 3:
            for transcript_id in eles[1].split(','):
                transcript_id_to_code[transcript_id] = eles[0]
    code_to_transcript_id_set = {'=': set()}
    for transcript_id in transcript_id_set:
        if transcript_id in transcript_id_to_code:
            code = transcript_id_to_code[transcript_id]
            if code in code_to_transcript_id_set:
                code_to_transcript_id_set[code].add(transcript_id)
            else:
                code_to_transcript_id_set[code] = {transcript_id}
        else:
            code_to_transcript_id_set['='].add(transcript_id)
    with open(new_code_num_txt, 'w') as fw:
        for code, id_set in code_to_transcript_id_set.items():
            fw.write('{}\t{}\t{}\n'.format(code, ','.join(id_set), len(id_set)))
