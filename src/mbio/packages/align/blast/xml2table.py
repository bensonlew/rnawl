# -*- coding: utf-8 -*-
# __author__ = 'sheng.he'

import os
import re
from Bio.Blast import NCBIXML
import pandas as pd
import numpy as np


default_header = ['Score', 'E-Value', 'HSP-Len', 'Identity-%', 'Similarity-%', 'Query-Name', 'Q-Len', 'Q-Begin',
                  'Q-End', 'Q-Frame', 'Hit-Name', 'Hit-Len', 'Hsp-Begin', 'Hsp-End', 'Hsp-Frame', 'Hit-Description']

all_values = ['Score', 'E-Value', 'HSP-Len', 'Identity-%', 'Similarity-%', 'Identity', 'Positives',
              'Query-Name', 'Q-Len', 'Q-Begin', 'Q-End', 'Q-Frame', 'Hit-Name', 'Hit-Len', 'Hsp-Begin',
              'Hsp-End', 'Hsp-Frame', 'Hit-Description', 'Q-Strand', 'Hsp-Strand', 'Mismatch', 'Gapopen_num']

# blast outfmt 6 默认格式，名字为blast中的默认名字
blast_detault = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
                 "send", "evalue", "bitscore"]
# 翻译blast输出默认名称为本模块默认的名称
transform_blast_default = ["Query-Name", "Hit-Name", "Identity-%", "HSP-Len", "Mismatch", "Gapopen_num",
                           "Q-Begin", "Q-End", "Hsp-Begin", "Hsp-End", 'E-Value', "Score"]


def xml2table(xml_fp, table_out, header=None, anno_head=True, hit_id=None):
    """

    :param xml_fp: 输入xml文件路径
    :param table_out: 输出文件路径
    :param header: 选择列，列的选择必须是上面定义的 all_values中的值组成的列表
    :param anno_head: 是否写入表头
    :param hit_id: 这里为兼容新版本的nt库的代码用，表示含义为是hit_id的取的方式改为从des中取
    """
    if header:
        for i in header:
            if i not in all_values:
                raise Exception('无法获取的值:{}\n可用的值:{}'.format(i, '\t'.join(all_values)))
    else:
        header = default_header
    if not os.path.isfile(xml_fp):
        raise Exception('输入xml文件不存在:{}'.format(xml_fp))
    with open(xml_fp) as f, open(table_out, 'w') as w:
        if anno_head:
            w.write('\t'.join(header) + '\n')
        records = NCBIXML.parse(f)
        values = {i: 'N/A' for i in all_values}
        for rec in records:
            query = re.split(' ', rec.query, maxsplit=1)[0]
            for align in rec.alignments:
                for hsp in align.hsps:
                    one_hsp = values.copy()
                    one_hsp['Query-Name'] = query
                    if hit_id:
                        one_hsp['Hit-Name'] = (align.hit_def).split(" ")[0]
                        one_hsp['Hit-Description'] = " ".join(align.hit_def.split(" ")[1:])
                    else:
                        one_hsp['Hit-Name'] = align.hit_id
                        one_hsp['Hit-Description'] = align.hit_def
                    one_hsp['Score'] = str(hsp.score)
                    one_hsp['E-Value'] = str(hsp.expect)
                    one_hsp['HSP-Len'] = str(hsp.align_length)
                    one_hsp['Identity'] = str(hsp.identities)
                    one_hsp['Positives'] = str(hsp.positives)
                    one_hsp['Q-Len'] = str(rec.query_length)
                    one_hsp['Q-Begin'] = str(hsp.query_start)
                    one_hsp['Q-End'] = str(hsp.query_end)
                    one_hsp['Q-Frame'] = str(hsp.frame[0])
                    one_hsp['Hit-Len'] = str(align.length)
                    one_hsp['Hsp-Begin'] = str(hsp.sbjct_start)
                    one_hsp['Hsp-End'] = str(hsp.sbjct_end)
                    one_hsp['Hsp-Frame'] = str(hsp.frame[1])
                    one_hsp['Q-Strand'] = str(hsp.strand[0])
                    one_hsp['Hsp-Strand'] = str(hsp.strand[1])
                    one_hsp['Identity-%'] = str(round(float(hsp.identities) / hsp.align_length, 3) * 100)
                    one_hsp['Similarity-%'] = str(round(float(hsp.positives) / hsp.align_length, 3) * 100)
                    one_hsp['Mismatch'] = str(int(hsp.align_length) - int(len(hsp.match)))  # 此处使用len(hsp.match) 还是hsp.indentities没有具体证据
                    one_hsp['Gapopen_num'] = str(hsp.gaps)  # 此处的gaps不是 gapopen 而是 gaps总数，因为xml获取不到
                    line = list()
                    for i in header:
                        line.append(one_hsp[i])
                    w.write('\t'.join(line) + '\n')
                    # line = ""
                    # for i in header:
                    #     line += one_hsp[i] + '\t'
                    # w.write(line.strip() + '\n')
    return table_out

def extract_xml(xml_file,output):
    """
    :param xml_file: 每个基因注释多条记录的xml文件
    :param output: 每个基因注释一条记录的xml文件
    :return:
    """
    ma =0
    with open (xml_file,'r') as f,open (output,'w') as p:
        files =f.readlines()
        for line in files:
            if re.search(r'\<\/Hit\>',line) and ma == 0:
                ma =1
                p.write(line)
            if re.search(r'\<\/Iteration_hits\>',line):
                ma =0
            if ma ==0:
                p.write(line)
            elif ma ==1:
                pass
    return output

def table_to_table(blast_table,output):
    """
    :param blast_table: blast outfmt 6表格格式的默认输出格式
    :param output: blast outfmt 6表格格式的默认输出格式
    :return: 主要检查outfmt 6表格格式第10列是否含gi号，如果有进行转化只有assesion号
    """
    with open(blast_table, 'r') as f, open(output, 'w') as p:
        files = f.readlines()
        p.write(files[0])
        for line in files[1:]:
            lines = line.split('\t')
            m = re.search('gi\|(.*)\|(.*)\|(.*)\|', lines[10])
            lines[10]= m.group(3)
            p.write('\t'.join(lines))
        return output

def xml2blast6(xml_fp, table_out, anno_head=False):
    """
    转换xml为blast outfmt 6表格格式的默认输出格式
    """
    xml2table(xml_fp, table_out, header=transform_blast_default, anno_head=anno_head)

def xml2table_coverage(xml_fp, table_out,head=None, anno_head=False,query=True,hit_ref=True):
    """
    :修改表头，计算覆盖度  zouguanqing 201903
    """
    if not head:
        head = ["Query-Name", "Hit-Name", "Identity-%", "HSP-Len", "Mismatch", "Gapopen_num",
               "Q-Begin", "Q-End", "Hsp-Begin", "Hsp-End", 'E-Value', "Score", "Q-Len","Hit-Len" ]
        if hit_ref == True:
            head[1] = "Hit-Description"
    tmp_table = table_out + '_tmp'
    xml2table(xml_fp, tmp_table, header=head, anno_head=True)

    f = pd.read_table(tmp_table,sep='\t')
    if query:
        f['coverage'] = np.abs((f['Q-End']-f['Q-Begin'])/f['Q-Len'])
    else:
        f['coverage'] = f['HSP-Len']/f['Hit-Len']
    if anno_head:
        heads = head.append('coverage')
        f.to_csv(table_out,sep='\t', index=False,header=heads)
    else:
        f.to_csv(table_out,sep='\t',index=False, header=False)




if __name__ == '__main__':  # for test
    # xml2table('C:\\Users\\sheng.he.MAJORBIO\\Desktop\\annotation\\annotation\\NR\\transcript.fa_vs_nr.xml',
            #   'C:\\Users\\sheng.he.MAJORBIO\\Desktop\\blast_test.xls', ['Score', 'E-Value', 'HSP-Len', 'Positives'])
    a = 'C:\\Users\\sheng.he.MAJORBIO\\Desktop\\fasta_5_vs_animal.xml'
    b = 'C:\\Users\\sheng.he.MAJORBIO\\Desktop\\blast_test_3.xls'
    xml2blast6(a, b)
