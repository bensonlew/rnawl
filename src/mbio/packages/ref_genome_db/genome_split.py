# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from BCBio import GFF
import sys
import os
import re
from optparse import OptionParser
from Bio import SeqIO


parser = OptionParser(description='split fasta too long')
parser.add_option('-g', '--gtf', dest='gtf', help='gtf file')
parser.add_option('-f', '--fasta', dest='fasta', help='fasta file')
parser.add_option('-n', '--nlength', dest='min_n_len', help='N length of split position', type=int, default = 10000)
parser.add_option('-m', '--max_length', dest='max_len', help='fasta maxlength', type=int, default = 256000000)
(opts, args) = parser.parse_args()

def uniq_gene(pos_list):
    # 去除重叠区域

    uniq_list = list()
    uniq_pos = pos_list[0]
    for pos in pos_list[1:]:
        if pos[0] <= uniq_pos[1]:
            if pos[1] > uniq_pos[1]:
                uniq_pos[1] = pos[1]
        else:
            uniq_list.append(uniq_pos)
            uniq_pos = pos
    uniq_list.append(uniq_pos)
    return uniq_list


def main(gtf_file, fasta_file, max_len=256000000, min_n_len=10000):
    print max_len, min_n_len
    gene_pos, chrom_set = get_gene_pos(gtf_file)
    # print gene_pos, chrom_set
    inter_pos = get_inter_pos(gene_pos, chrom_set)

    # print inter_pos

    split(inter_pos, fasta_file, max_len, min_n_len, gtf=gtf_file)


def get_n_dict(seq):
    n_pos_dict = dict()
    range_n = False
    pos_idx = 1
    for i in range(len(seq)):
        if seq[i].lower() == "n":
            if not range_n:
                n_pos_dict[pos_idx] = [i, 1]
            else:
                n_pos_dict[pos_idx][1] += 1
            range_n = True
        else:
            if range_n:
                pos_idx += 1
            range_n = False
    return n_pos_dict

def get_n_pos(min_n_len, intron_seq, start, n_path=None):
    # 查找连续N的位置并返回中间点的坐标
    if len(intron_seq) < min_n_len:
        return False
    else:
        n_pos_dict = get_n_dict(intron_seq)
        if len(n_pos_dict) == 0:
            return False
        largest_n = sorted(n_pos_dict.values(), key=lambda x:x[1], reverse=True)[0]
        if largest_n[1] > min_n_len:
            return start + largest_n[0] + largest_n[1]/2
        else:
            return False
        '''
        # n_path = re.compile(".*(" + "N" * min_n_len + ").*", re.I)
        # print intron_seq
        # print n_path
        # print "type is {}".format(type(intron_seq))
        m = re.search(n_path, intron_seq)
        if m:
            return start + m.start(1) + min_n_len/2
        else:
            return False
        '''


def seq_split(seq, pos, start, max_len, min_n_len, break_list, gene_dis):
    # 分割单条序列
    break_pre = 0
    for i,apos in enumerate(pos):
        if apos[1] > max_len + start:
            break_pre = i-1
            break

    # print "pos is start {}".format(start)
    # n_path = re.compile(".*(" + "N" * min_n_len + ").*", re.I)
    break_point = False

    # print "i is {}".format(i)
    for j in range(i-1, 0, -1):
        if pos[j][1] - gene_dis - (pos[j][0] + gene_dis) < min_n_len:
            continue
        intron_seq = seq[pos[j][0] - 1 + gene_dis: pos[j][1] - gene_dis]
        # print intron_seq
        start = pos[j][0] - 1 + gene_dis
        break_point = get_n_pos(min_n_len, intron_seq, start)
        if seq[pos[j][0] - 1] < start:
            break
        if break_point:
            break
        else:
            pass

    if not break_point:
        raise Exception("无法找到合适的切割位点")
    else:
        break_list.append(break_point)

    if len(seq) - break_point < max_len:
        return break_list
    else:
        seq_split(seq, pos, break_point, max_len, min_n_len, break_list, gene_dis)

    return break_list


def split(inter_pos, fasta, max_len, min_n_len, gene_dis=500, gtf=None):
    # 分割序列
    pos_dict = dict()
    with open(fasta + "split.fa", 'w') as fo, open(fasta + "split.txt", 'w') as fjo:
        pos_json = dict()
        for seq in SeqIO.parse(fasta, "fasta"):
            seq_str = seq.seq
            if str(seq.id) not in inter_pos:
                continue
            pos = inter_pos[str(seq.id)]
            break_list = list()
            print len(seq_str)
            if len(seq_str) > int(max_len):
                break_list = seq_split(str(seq_str), pos, 0, max_len, min_n_len, break_list, gene_dis)
                # print seq
                print "break list is {}".format(break_list)
                start = 0
                for break1 in break_list:
                    fo.write(">{}\n{}\n".format(str(seq.id) + "__" + str(start),
                                                str(seq.seq)[start:break1]))
                    fjo.write("\t".join([str(seq.id) + "__" + str(start), str(seq.id), str(start + 1), str(break1)]) + "\n")
                    if seq.id not in pos_dict:
                        pos_dict[seq.id] = [[str(seq.id) + "__" + str(start), start, break1]]
                    else:
                        pos_dict[seq.id].append([str(seq.id) + "__" + str(start), start, break1])
                    start = break1

                fo.write(">{}\n{}\n".format(str(seq.id) + "__" + str(start),
                                                str(seq.seq)[start:len(seq_str)]))
                fjo.write("\t".join([str(seq.id) + "__" + str(start), str(seq.id), str(start + 1 ), str(len(seq_str))]) + "\n")
                if seq.id not in pos_dict:
                    pos_dict[seq.id] = [[str(seq.id) + "__" + str(start), start, len(seq_str)]]
                else:
                    pos_dict[seq.id].append([str(seq.id) + "__" + str(start), start, len(seq_str)])

            else:
                fo.write(">{}\n{}\n".format(str(seq.id), str(seq.seq)))
                fjo.write("\t".join([str(seq.id), str(seq.id), str(1), str(len(seq_str))]) + "\n")
                pos_dict[seq.id] = [[str(seq.id), 0, len(seq_str)]]

        print pos_dict
        with open(gtf, 'r') as f_in, open(gtf + "split.txt", 'w') as fgo:
            for line in f_in:
                cols = line.split("\t")
                if len(cols) < 9:
                    fgo.write(line)
                else:
                    if cols[0] in pos_dict and len(pos_dict[cols[0]]) >= 2:
                        for sub_pos in pos_dict[cols[0]]:
                            if int(cols[3]) > sub_pos[1] and int(cols[4]) <= sub_pos[2]:
                                cols[0] = sub_pos[0]
                                cols[3] = str(int(cols[3]) - sub_pos[1])
                                cols[4] = str(int(cols[4]) - sub_pos[1])
                                fgo.write("\t".join(cols))
                                break
                    else:
                        fgo.write(line)



def get_inter_pos(gene_pos, chrom_set):
    # 获取每条染色体基因间区坐标
    inter_pos_dict = dict()
    for chrom in chrom_set:
        gene_regions = [[g['start'], g['end']]  for g in gene_pos.values() if g['chrom'] == chrom]
        gene_regions_sorted = sorted(gene_regions, key=lambda x:x[0])
        gene_regions_uniq = list()
        uniq_list = uniq_gene(gene_regions_sorted)
        inter_pos = list()
        for i in range(1, len(uniq_list)):
            inter_pos.append([uniq_list[i-1][1], uniq_list[i][0]])
        inter_pos_dict[chrom] = inter_pos

    return inter_pos_dict


def get_gene_pos(gtf_file):
    # 获取基因所在位置信息
    gene_pos = dict()
    chrom_set = set()
    pattern1 = re.compile("transcript_id \"(.+)\";.*gene_id \"(.+)\"", re.I)
    pattern2 = re.compile("gene_id \"(.+)\";.*transcript_id \"(.+)\"", re.I)
    with open(gtf_file, 'r') as in_handle:
        for line in in_handle.readlines():

            cols = line.strip("\n").split("\t")
            # print cols[-1]
            tran_id = ""
            gene_name = ""

            m = re.match(pattern1, cols[-1])
            if m:
                tran_id = m.group(1)
                gene_id = m.group(2)
                start = int(cols[3])
                end = int(cols[4])

                if gene_id in gene_pos:
                    if gene_pos[gene_id]['start'] > start:
                        gene_pos[gene_id]['start'] = start
                    if gene_pos[gene_id]['end'] < end:
                        gene_pos[gene_id]['end'] = end
                else:
                    gene_pos[gene_id] = {
                        'start': start,
                        'end': end,
                        'chrom': cols[0]
                    }
                    chrom_set.add(cols[0])

            else:
                m = re.match(pattern2, cols[-1])
                if m:
                    tran_id = m.group(2)
                    gene_id = m.group(1)
                    start = int(cols[3])
                    end = int(cols[4])
                    if gene_id in gene_pos:
                        if gene_pos[gene_id]['start'] > start:
                            gene_pos[gene_id]['start'] = start
                        if gene_pos[gene_id]['end'] < end:
                            gene_pos[gene_id]['end'] = end
                    else:
                        gene_pos[gene_id] = {
                            'start': start,
                            'end': end,
                            'chrom': cols[0]
                        }
                    chrom_set.add(cols[0])
    return gene_pos, chrom_set


if __name__ == '__main__':
    if opts.gtf and opts.fasta:
        main(opts.gtf, opts.fasta, max_len=opts.max_len, min_n_len=opts.min_n_len)
    else:
        parser.print_help()
