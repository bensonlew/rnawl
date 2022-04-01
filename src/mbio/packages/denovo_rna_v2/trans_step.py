# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

from Bio import SeqIO
import fileinput
import re
import os
import subprocess
import urllib2
from collections import defaultdict
import regex


def step_count(fasta_file, fasta_to_txt, group_num, step, stat_out, min_len=200):
    """
    步长统计
    :param fasta_file: 输入的fa文件
    :param fasta_to_txt:输出的统计数据txt
    :param group_num:按照步长统计几组
    :param step:步长
    :param stat_out:统计的数据汇总信息txt（fasta_to_txt文件的汇总）
    :return:
    """
    with open(fasta_to_txt, "w") as f:
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            ID = seq_record.description.strip().split(" ")[0]
            new_trans_list = ID + "\t" + str(len(seq_record)) + "\n"
            f.write(new_trans_list)
    with open(fasta_to_txt, "r") as r, open(stat_out, "a") as w:
        sample_name = os.path.basename(fasta_file).split('.fa')[0]
        w.write("Length" + "\tNumber" + "\tPercent"+ "\n")
        trans_list = []
        amount_group = []
        element_set = set("")
        for line in r:
            line = line.strip().split("\t")
            number = line[1]
            trans_list.append(number)
        for f in trans_list:
            for i in range(group_num):
                if (int(f) >= (i * step)) and (int(f) < ((i+1) * step)):
                    amount_group.append(i)
                else:
                    pass
                    # amount_group.append(group_num+1)
                element_set.add(i)
        amount_group.sort()
        top_sum = 0
        all_sum = sum([amount_group.count(i) for i in element_set])

        for i in element_set:
            num_statistics = amount_group.count(i)
            pct = format(float(num_statistics)/all_sum, '.2%')
            if str(i) == '0':
                area_line = str(min_len) + "~" + str((i + 1) * step) + "\t" + str(num_statistics) + "\t" + str(pct) + "\n"
                w.write(area_line)
                top_sum += int(num_statistics)
            elif i < (group_num-1):
                area_line = str(i * step + 1) + "~" + str((i+1) * step) + "\t" + str(num_statistics) + "\t" + str(pct) + "\n"
                w.write(area_line)
                top_sum += int(num_statistics)
            else:
                pct = format(float(str(len(trans_list)-int(top_sum)))/all_sum, '.2%')
                area_line = ">" + str(i * step) + "\t" + str(len(trans_list)-int(top_sum)) + "\t" + str(pct) + "\n"
                end_line = "total" + "\t" + str(len(trans_list)) + "\t" + '100%'+ "\n"
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
                        if candidateDic. has_key(tagpair[1]):
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
        group = max_num/step
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
            area_line = str(n * step) + "~" + str((n + 1) * step) + "\t" + str(final_value) + "\t" + str(final_ids) + "\n"
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
            if str(method) == "stringtie":
                if seq_type.split("\t")[0] == "exon":
                    trans_dic[txpt_id] += 1
            elif str(method) == "cufflinks":
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
if __name__ == '__main__':
#     merged_gtf = '/mnt/ilustre/users/sanger-dev/workspace/20170401/Single_assembly_module_tophat_cufflinks/Assembly/output/Cuffmerge/merged.gtf'
#     old_genes_gtf = "/mnt/ilustre/users/sanger-dev/workspace/20170407/Single_assembly_module_tophat_stringtie_gene2/Assembly/output/NewTranscripts/old_genes.gtf"
#     new_genes_gtf = "/mnt/ilustre/users/sanger-dev/workspace/20170407/Single_assembly_module_tophat_stringtie_gene2/Assembly/output/NewTranscripts/old_genes.gtf"
#     output1 = "/mnt/ilustre/users/sanger-dev/workspace/20170401/Single_assembly_module_tophat_cufflinks/Assembly/output/Cuffmerge/1.txt"
#     output2 = "/mnt/ilustre/users/sanger-dev/workspace/20170401/Single_assembly_module_tophat_stringtie_gene2/Assembly/output/assembly_newtranscripts/2.txt"
#     output3 = "/mnt/ilustre/users/sanger-dev/workspace/20170401/Single_assembly_module_tophat_stringtie_gene2/Assembly/output/assembly_newtranscripts/3.txt"
#     output4 = "/mnt/ilustre/users/sanger-dev/workspace/20170401/Single_assembly_module_tophat_stringtie_gene2/Assembly/output/assembly_newtranscripts/4.txt"
#     class_code_count(merged_gtf, output1)
#     merged_gtf = "O:\\Users\\zhaoyue.wang\\Desktop\\merged1.gtf"
#     output1 = "/mnt/ilustre/users/sanger-dev/workspace/20170410/Single_assembly_module_tophat_stringtie_zebra/Assembly/output/Statistics/1.txt"
#     output2 = "/mnt/ilustre/users/sanger-dev/workspace/20170410/Single_assembly_module_tophat_stringtie_zebra/Assembly/output/Statistics/2.txt"
#     output3 = "O:\\Users\\zhaoyue.wang\\Desktop\\3.txt"
#     output4 = "O:\\Users\\zhaoyue.wang\\Desktop\\4.txt"
#     gene_trans_exon(merged_gtf, "stringtie", output3, output4)
    # gene_trans_exon(merged_gtf, output1, output2)
    # gtf_file = "O:\\Users\\zhaoyue.wang\\Desktop\\old_trans.gtf.trans"
    # gtf_files = "/mnt/ilustre/users/sanger-dev/workspace/20170410/Single_assembly_module_tophat_stringtie_zebra/Assembly/assembly_newtranscripts/old_trans.gtf.trans"
    # count_trans_or_exons(gtf_files, 5, output1, output2)
    fa_file = '/mnt/ilustre/users/sanger-dev/sg-users/shijin/Refrna_demo1/RefrnaAssemble/output/StringtieMerge/change_id_merged.fa'
    txt_file = '/mnt/ilustre/users/sanger-dev/sg-users/wangzhaoyue/moduletest/200.txt'
    step_count(fa_file, txt_file, 10, 200, "/mnt/ilustre/users/sanger-dev/sg-users/wangzhaoyue/moduletest/final_200.txt")
