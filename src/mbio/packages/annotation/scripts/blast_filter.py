# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2017.06.09
import sys


def blast_filter(blastout_table, blast_fliter, evalue, score, similarity, identity):
    """进行blast结果筛选，参数：evalue, score, similarity, identity"""
    evalue_list, score_list, similarity_list, identity_list = [], [], [], []
    with open(blastout_table, "rb") as f, open(blast_fliter, "wb") as w:
        lines = f.readlines()
        w.write(lines[0])
        for i in range(1, len(lines)):
            line = lines[i].strip().split("\t")
            f_score = line[0]
            f_evalue = line[1]
            f_identity = line[3]
            f_similarity = line[4]
            if evalue >= float(f_evalue):
                evalue_list.append(i)
            if score <= float(f_score):
                score_list.append(i)
            if similarity <= float(f_similarity):
                similarity_list.append(i)
            if identity <= float(f_identity):
                identity_list.append(i)
        for j in range(1, len(lines)):
            if j in evalue_list:
                if j in score_list:
                    if j in similarity_list:
                        if j in identity_list:
                            w.write(lines[j])


blast_filter(blastout_table=sys.argv[1], blast_fliter=sys.argv[2], evalue=float(sys.argv[3]), score=float(sys.argv[4]), similarity=float(sys.argv[5]), identity=float(sys.argv[6]))
