#!/usr/bin/python
# -*- coding: utf-8 -*-
import re
import os

def blastStat(blast_infile, blast_outfile, Score = None, Evalue = None, HSP_len = None,Identity = None, Similarity = None):
    """
    对blast结果E-value、Score、HSP-Length、Identity、Similarity进行统计  16列
    :param: blast_infile, blast输入文件 xls格式
    :param: blast_outfile, blast输出文件 xls格式
    :param: Score, Evalue , HSP_len ,Identity , Similarity 为浮点型或整型区间
    """
    if not os.path.exists(blast_infile):
        raise Exception("{}文件不存在！".format(blast_infile))
    query = {}
    with open(blast_infile, 'r+') as f1:
        names=f1.readline().strip().split("\t")
        print names
        for lines in f1:
            print lines
            line=lines.strip().split("\t")
            query_id = line[5]
            if query_id not in query.keys():
                query[query_id] = {}
                query[query_id]["score"] = [line[0]]
                query[query_id]["evalue"] = [line[1]]
                query[query_id]["hspLen"] = [line[2]]
                query[query_id]["identity"] = [line[3]]
                query[query_id]["similarity"] = [line[4]]
                query[query_id]["query"] = [line]
            else:
                query[query_id]["score"].append(line[0])
                query[query_id]["evalue"].append(line[1])
                query[query_id]["hspLen"].append(line[2])
                query[query_id]["identity"].append(line[3])
                query[query_id]["similarity"].append(line[4])
                query[query_id]["query"].append(line)
    print query
    with open(blast_outfile, 'w+') as f2:
        for keys in query.keys():
            data = query[keys]
            min_value = min(data["evalue"])
            minEvalue = data["evalue"][data["evalue"].index(min_value)]
            print minEvalue, Evalue[1]
            print Evalue[0]
            print minEvalue > Evalue[1]
            if Evalue:
                    """
                    if not (minEvalue < Evalue[1] and minEvalue > Evalue[0]):
                        print "haha"
                        print "{}没有找到最佳序列".format(keys)
                        pass
                    else:
                    """
                    print 'heihei'
                    count = 0
                    count_data = []
                    #统计最小evalue值数目
                    for t in range(len(query[keys]["evalue"])):
                        if query[keys]["evalue"][t] == min_value:
                            count += 1
                            count_data.append(t)
                    #如果最小evalue值有至少两个，则比较score大小，取score较大者
                    print count
                    print count_data
                    if count > 1:
                        new_data = {}
                        new_data['score'] = []
                        new_data['query'] = []
                        new_data['hspLen'] = []
                        new_data['identity'] = []
                        new_data['similarity'] = []
                        for ss in count_data:
                            new_data["score"].append(data["score"][ss])
                            new_data["query"].append(data["query"][ss])
                            new_data["hspLen"].append(data["hspLen"][ss])
                            new_data["identity"].append(data["identity"][ss])
                            new_data["similarity"].append(data["similarity"][ss])
                        print new_data
                        print "evalue值有相等的！"
                        #默认当evalue存在相等的,但score不会相等
                        score = max(new_data["score"])
                        score_index = new_data["score"].index(score)
                        _line = new_data['query'][score_index]
                        hsp = new_data["hspLen"][score_index]
                        identity = new_data["identity"][score_index]
                        similarity = new_data["similarity"][score_index]
                    else:
                        _line = data['query'][data['evalue'].index(min_value)]
                        score = data['score'][data['evalue'].index(min_value)]
                        hsp = data['hspLen'][data['evalue'].index(min_value)]
                        identity = data['identity'][data['evalue'].index(min_value)]
                        similarity = data['similarity'][data['evalue'].index(min_value)]
                    #依次判断Score, HSP_len, Identity, Similarity 是否满足设定的条件
                    if Score:
                        if score_max < Score[0] or score_max > Score[1]:
                            next
                        else:
                            pass
                    if HSP_len:
                        if hsp > HSP_len[1] or hsp < HSP_len[0]:
                            next
                        else:
                            pass
                    if Identity:
                        if identity > Identity[1] or identity < Identity[0]:
                            next
                        else:
                            pass
                    if Similarity:
                        if similarity > Similarity[1] or similarity < Similarity[0]:
                            next
                        else:
                            pass
                    f2.write("\t".join(_line)+"\n")
    print 'end'

if __name__ == "__main__":
    #blastStat("blast.xls","blast_new.xls",Evalue = [1e-100,1e-10])
    blastStat("/mnt/ilustre/users/sanger-dev/sg-users/konghualei/denovo_rna/annotation/blast.xls", "blast_new.xls", Evalue = [1e-100,1e-10])