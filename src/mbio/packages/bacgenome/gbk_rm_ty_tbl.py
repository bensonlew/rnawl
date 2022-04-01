#!/sur/bin/python
#-*- coding: utf-8 -*-
#Deal gbk file
#20181225

import sys
import re
import copy

def deal_tbl(tbl):
    f = open(tbl, 'r')
    seqs = dict()
    name = ''
    loc = ''
    tmp = {}
    mk_cds = 0
    mk_rna = 0
    lines = f.readlines()
    gene_pat = re.compile('\d+\t\d+\tgene')
    cds_pat = re.compile('\d+\t\d+\tCDS')
    rrna_pat = re.compile('\d+\t\d+\trRNA')
    trna_pat = re.compile('\d+\t\d+\ttRNA')
    info_pat = re.compile('^([^\t]*)\t([^\t]*)$')
    for i in lines:
        i = i.strip()        
        m = re.match('>Feature\s*(.*)', i)
        if m:
            name = m.group(1)
            seqs[name] = dict()
        if re.match(gene_pat, i):
            if loc != '':
                seqs[name][loc] = tmp
            mk_cds = 0
            tmp = {}
        if re.match(cds_pat, i):
            mk_cds = 1
        if re.match(rrna_pat, i) or re.match(trna_pat, i):
            mk_rna =1
        m = re.match(info_pat, i)
        if m and len(seqs) == 1:
            #tmp[m.group(1)] = m.group(2)
            tmp[m.group(1)] = m.group(2).strip() #zouguanqing 20190716
            if m.group(1) == 'locus_tag':
                loc = m.group(2)
        if m and mk_cds == 1:
            #tmp[m.group(1)] = m.group(2)
            tmp[m.group(1)] = m.group(2).strip() #zouguanqing 20190716
            if m.group(1) == 'locus_tag':
                loc = m.group(2)
        if m and mk_rna == 1:
            tmp[m.group(1)] = m.group(2).strip()  #zouguanqing 20190716
            if m.group(1) == 'locus_tag':
                loc = m.group(2)
    seqs[name][loc] = tmp

    return seqs

def deal_gbk(gbk):
    f = open(gbk, 'r')
    lines = f.readlines()
    scaf_pat = re.compile('^LOCUS\s*(\S+)')
    gene_pat = re.compile('^\s+gene\s+\S*$')
    cds_pat = re.compile('^\s+CDS\s+\S*$')
    rrna_pat = re.compile('^\s+rRNA\s+\S*$')
    trna_pat = re.compile('^\s+tRNA\s+\S*$')
    source_pat = re.compile('^\s+source\s+\S*$')
    loc_pat = re.compile('\s*/locus_tag="(.*)"')
    condon_pat = re.compile('\s*/codon_start=1\s+S*$')
    tra_pat = re.compile('^\s*/translation=".*\s+')
    lef_pat = re.compile('^\s[A-Z]*"\s+')
    pro_pat = re.compile('\s*/product="(.*)"')
    start = ''
    end = ''
    mid = list()
    s = 1
    e = 0
    tmp = ''
    source_mk = 0
    gene_mk = 0
    cds_mk = 0
    rrna_mk = 0
    ret_dic = dict()
    name_list = []
    name = ''
    tra_pro =''
    ss = 0
    gene_s = ''
    tmp_list = []
    for i in lines:  #遍历gbk
        m_scaf = re.match(scaf_pat,i)
        if m_scaf:
            if name != '':
                ret_dic[name] = [start, end, mid]
                name_list.append(name)
            name = m_scaf.group(1)
        if re.match('FEATURES', i):
            s = 0
        elif re.match('ORIGIN', i):
            e = 1
        if s == 1:
            start += i
        elif e == 1 :
            end += i
        elif s == 0  and e == 0:
            if re.match(source_pat, i):
                source_mk = 1
                gene_mk = 0
                cds_mk = 0
                tmp = ''
                tmp_list.append(tmp)

            if re.match(gene_pat, i):
                tmp_list.append(tmp)
                tmp_list.insert(0, gene_s)
                mid.append(tmp_list)
                gene_mk = 1
                source_mk = 0
                cds_mk = 0
                rrna_mk = 0
                tmp = {'content':'','loc':'','product':''}
                gene_s = ''
                tmp_list = []

            if gene_mk == 0:
                if re.match(cds_pat, i):
                    continue
                elif re.match(pro_pat, i):
                    continue
                elif re.match(condon_pat, i):
                    continue
                if re.match(tra_pat, i):
                    ss =1
                    continue
                elif ss ==1 and re.match(lef_pat,i):
                    continue
                    ss = 0
                elif ss ==1:
                    continue
            elif gene_mk ==1:
                if re.match(cds_pat, i) or re.match(trna_pat, i) or re.match(rrna_pat, i):
                    if tmp['content']!='':
                        tmp_list.append(tmp)
                    tmp = {'content':'','loc':'','product':''}
                    cds_mk = 1

            if source_mk :
                tmp += i
            elif gene_mk ==0:
                if cds_mk == 1:
                    tmp['content'] += i
                    m_loc =  re.match(loc_pat, i)
                    m_pro = re.match(pro_pat, i)
                    if m_loc:
                        tmp['loc'] = m_loc.group(1)
                    if m_pro:
                        tmp['product'] = m_pro.group(1).strip()  #zouguanqing 20190716
                else:
                    gene_s += i
            elif gene_mk ==1:
                if cds_mk == 1:
                    tmp['content'] += i
                    m_loc =  re.match(loc_pat, i)
                    m_pro = re.match(pro_pat, i)
                    if m_loc:
                        tmp['loc'] = m_loc.group(1)
                    if m_pro:
                        tmp['product'] = m_pro.group(1).strip()  #zouguanqing 20190716
                else:
                    gene_s += i
    tmp_list.append(tmp)
    tmp_list.insert(0, gene_s)
    mid.append(tmp_list)
    ret_dic[name] = [start,end,mid]
    name_list.append(name)

    return ret_dic, name_list

def deal_middle_part(mid, search):
    #new_mid = copy.deepcopy(mid)
    ret_s = ''
    for i in mid:
        #print i
        if i == '':
            continue
            #new_mid.remove(i)

        elif isinstance(i, list):
            dic =[]  # zouguanqing 20190716
            for j in i:
                #dic =[]
                if isinstance(j, dict):
                    dic.append(j)
            if len(dic) ==1:
                for j in i:
                    if isinstance(j, dict):
                        ret_s += j['content']
                    else:
                        ret_s += j
            elif len(dic) >1:
                mark_do = 0  #limit to  produce one cds
                for j in i:
                    if isinstance(j, dict):
                        loc = j['loc']

                        if loc in search.keys():
                            if j['product'] == search[loc]['product'] and mark_do == 0:
                                ret_s += j['content']
                                mark_do = 1
                    else:
                        ret_s += j

            elif len(dic) == 0:
                for j in i:
                    ret_s += j

        else:
            ret_s +=i
    return ret_s

def produce_gbk(ret_gbk,name_list, ret_tbl, out):
    fw = open(out, 'w')
    for k in name_list:
        if k not in ret_tbl.keys():
            raise Exception('KEY ERROR')
        mid = ret_gbk[k][2]
        tbl = ret_tbl[k]
        ret_str = deal_middle_part(mid, tbl)

        gbk = ret_gbk[k][0] + ret_str + ret_gbk[k][1]
        fw.write(gbk)


if __name__ == '__main__':
    tbl = sys.argv[1]
    gbk = sys.argv[2]
    out = sys.argv[3]
    ret_tbl = deal_tbl(tbl)
    ret_gbk, name_list = deal_gbk(gbk)
    produce_gbk(ret_gbk, name_list, ret_tbl,out)














    
