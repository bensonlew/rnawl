# -*- coding:utf-8 -*-
# !/usr/bin/python
# zouguanqing
# 20190218

import sys
import copy
import argparse


def type_sort():
    num1 = 'glimmer'
    num2 = 'genemarks'
    num3 = 'prodigal'
    return [num1, num2, num3]


def get_type_info(infile, name_index, pos1_index, pos2_index, type, strand_pos=None, header=0, seq_index=None,
                  seq_split=False):
    ret = {}
    if seq_index == None:
        ret['-'] = []
    with open(infile) as fr:
        while header > 0:
            fr.readline()
            header -= 1
        for line in fr:
            line = line.strip()
            if line == '':
                continue
            sp = line.split('\t')
            name = type + ':' + sp[name_index]
            pos1 = int(sp[pos1_index])
            pos2 = int(sp[pos2_index])
            if pos1 > pos2:
                tmp = [name, pos2, pos1]
            else:
                tmp = [name, pos1, pos2]
            if strand_pos:
                tmp.append(sp[strand_pos])
            else:
                tmp.append(' ')
            if seq_index != None:
                if seq_split:
                    k = sp[seq_index].split(' ')[0]
                else:
                    k = sp[seq_index]

                if k not in ret.keys():
                    ret[k] = [tmp]
                else:
                    ret[k].append(tmp)
            else:
                ret['-'].append(tmp)
    return ret


def add_element_each(core, out, overlap, is_sort=None, rm_ori=False):
    if not rm_ori:
        new_core = copy.deepcopy(core)
    else:
        new_core = []

    if not is_sort:
        core = sorted(core, key=lambda a: a[1], reverse=False)
        out = sorted(out, key=lambda a: a[1], reverse=False)

    core_len = len(core)
    out_len = len(out)
    new_start = 0
    for e in out:
        p1 = e[1]
        p2 = e[2]
        start = new_start
        for i in range(start, core_len):
            current = core[i]
            c_p1 = current[1]
            c_p2 = current[2]
            if p2 < c_p1:
                new_core.append(e)
                new_start = i
                break
            elif p1 > c_p2:
                if i == core_len - 1:
                    new_core.append(e)
                continue
            elif p1 >= c_p1 and p2 <= c_p2:
                new_start = i
                break
            elif p2 - c_p1 + 1 >= overlap and p1 < c_p1:
                new_start = i
                break
            elif p2 - c_p1 + 1 < overlap and p1 < c_p1:
                new_core.append(e)
                new_start = i
                break
            elif c_p2 - p1 + 1 >= overlap and p2 > c_p2:
                new_start = i
                break
            elif c_p2 - p1 + 1 < overlap and p2 > c_p2:
                new_core.append(e)
                new_start = i
                break
            else:
                print "other condition: "
                print e
                print current
                exit()
    print('file1_gene_num: ' + str(core_len))
    print('file2_gene_num: ' + str(out_len))
    print('new_gene_num: ' + str(len(new_core)))
    return new_core


def add_element(core_dic, out_dic, overlap, is_sort=None, rm_ori=False):
    new_dic = {}
    core_k = core_dic.keys()
    out_k = out_dic.keys()
    core_k.extend(out_k)
    k_set = set(core_k)

    for k in sorted(k_set):
        print('result:' + k)
        if k not in core_dic.keys():
            new_dic[k] = out_dic[k]
            print('file1 has not ' + k + ' result.so add file2 ' + k + ' result, {} genes'.format(len(out_dic[k])))

        elif k not in out_dic.keys():
            print('file2 has not ' + k + ' result')
            if not rm_ori:
                new_dic[k] = core_dic[k]
        else:
            each_new = add_element_each(core_dic[k], out_dic[k], overlap, is_sort, rm_ori)
            new_dic[k] = each_new

    return new_dic


def main_(type_file_list, trna=None, rrna=None):
    all_type_sort = type_sort()
    current_types_sort = []
    for i in all_type_sort:
        for j in type_file_list:
            if i == j[0]:
                current_types_sort.append(j)

    c_index_info = current_types_sort[0][2]
    type = current_types_sort[0][0]
    infile = current_types_sort[0][1]

    gene_core = get_type_info(infile, c_index_info['name'], c_index_info['pos1'], c_index_info['pos2'],
                              type, c_index_info['strand'], header=1, seq_index=c_index_info['seq_id'], seq_split=True)
    print('####The starting file, {} : {} '.format(type, infile))
    if len(current_types_sort) > 2:
        for t, f, i in current_types_sort[1:]:
            print('####start merge {}: {}'.format(t, f))
            if t == 'prodigal':
                split_or = False
            else:
                split_or = True
            other = get_type_info(f, i['name'], i['pos1'], i['pos2'], t, i['strand'], header=1, seq_index=i['seq_id'],
                                  seq_split=split_or)
            gene_core = add_element(gene_core, other, 50)

    if trna or rrna:
        if trna:
            trna_dic = get_type_info(trna[1], trna[2]['name'], trna[2]['pos1'], trna[2]['pos2'], trna[0],
                                     header=1, seq_index=trna[2]['seq_id'], seq_split=True)

            print('####start merge {}: {}'.format('trna', trna[1]))
            gene_core = add_element(trna_dic, gene_core, 1, rm_ori=True)
        if rrna:
            rrna_dic = get_type_info(rrna[1], rrna[2]['name'], rrna[2]['pos1'], rrna[2]['pos2'], rrna[0],
                                     header=1, seq_index=rrna[2]['seq_id'], seq_split=True)

            print('####start merge {}: {}'.format('rrna', rrna[1]))
            gene_core = add_element(rrna_dic, gene_core, 1, rm_ori=True)

    return gene_core


def write_out(gene_data, prefix_gene, out_file='pos.predict', num_len=4):
    with open(out_file, 'w') as fw:
        for scaf in gene_data.keys():
            fw.write('>' + scaf + '\n')
            tmp = sorted(gene_data[scaf], key=lambda a: a[0])
            for id, ge in enumerate(tmp, 1):
                name = prefix_gene + '0' * (num_len - len(str(id))) + str(id)
                if ge[3] == '-' and ge[2] < ge[3]:
                    fw.write(name + '\t' + '\t'.join([str(ge[2]), str(ge[1]), ge[3]]) + '\n')
                else:
                    fw.write(name + '\t' + '\t'.join([str(ge[1]), str(ge[2]), ge[3]]) + '\n')


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('-p', '--prodigal', help='GFF file of prodigal')
    p.add_argument('-m', '--genemarks', help='GFF file of genemarks')
    p.add_argument('-g', '--glimmer', help='GFF file of glimmer')
    p.add_argument('-t', '--trna', help='GFF file of tRNA')
    p.add_argument('-r', '--rrna', help='GFF file of rRNA')
    a = p.parse_args()
    type_list = []
    index_info = {  ##指定各输入文件的信息所在的列。。name is not importance。
        'prodigal': {'name': 0, 'pos1': 3, 'pos2': 4, 'strand': 6, 'seq_id': 0},
        'glimmer': {'name': 0, 'pos1': 2, 'pos2': 3, 'strand': 4, 'seq_id': 1},
        'genemarks': {'name': 0, 'pos1': 2, 'pos2': 3, 'strand': 4, 'seq_id': 1},
        'trna': {'name': 0, 'pos1': 2, 'pos2': 3, 'seq_id': 1},
        'rrna': {'name': 0, 'pos1': 3, 'pos2': 4, 'seq_id': 1}
    }

    if a.prodigal:
        with open(a.prodigal) as f, open('prodigal_tmp.gff', 'w') as fw:
            for line in f:
                if line[0] == '#':
                    continue
                spline = line.split('\t')
                if spline[2] == 'CDS':
                    fw.write(line)
        type_list.append(['prodigal', 'prodigal_tmp.gff', index_info['prodigal']])
    if a.glimmer:
        type_list.append(['glimmer', a.glimmer, index_info['glimmer']])
    if a.genemarks:
        type_list.append(['genemarks', a.genemarks, index_info['genemarks']])

    trna = None
    rrna = None
    if a.trna:
        trna = ['trna', a.trna, index_info['trna']]
    if a.rrna:
        rrna = ['rrna', a.rrna, index_info['rrna']]
    gene_data = main_(type_list, trna, rrna)
    write_out(gene_data, 'orf')



