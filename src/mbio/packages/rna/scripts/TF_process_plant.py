# -*- coding: utf-8 -*-
# author: Moli
# last update: 20161025
#new database using hmmscan

#python TF_process_animal_v1.py ../hmmerdb/plant/planttfdb.hmm
# /mnt/ilustre/users/sanger-dev/app/database/refGenome/Animal/Primates/human/Homo_sapiens.GRCh38.pep.all.fa
# family_vs_DBD_plant.txt
# /mnt/ilustre/users/sanger-dev/workspace/20161116/Single_express_sample_v8_change_sample_name
# /Express/output/diff/featurecounts_count.txt.A_vs_B.edgeR.DE_results_name

import re
import sys
import os
import subprocess

'''
利用hmmer软件将未知序列比对到蛋白域中去
hmmscan <hmmdb> (seq)
'''
hmmer_path = '/mnt/ilustre/users/sanger-dev/app/bioinfo/align/hmmer-3.1b2-linux-intel-x86_64/binaries/'
# ref = '/mnt/ilustre/users/sanger-dev/sg-users/zhoumoli/TF/hmmerdb/planttfdb.hmm'
# query = '/mnt/ilustre/users/sanger-dev/sg-users/zhoumoli/TF/code/test_plant.fas'

cmd = "{}hmmscan {} {}".format(hmmer_path, sys.argv[1], sys.argv[2])

try:
    # subprocess.check_output(cmd, shell=True)
    res = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    mi = res.communicate()
    out = list(mi)
    out = out[0].split('\n')
except subprocess.CalledProcessError:
    raise Exception("运行出错！")

# family_db = '/mnt/ilustre/users/sanger-dev/sg-users/zhoumoli/TF/code/family_DBD.txt'
family_db = sys.argv[3]
query = []
phmmer = []
domain = []
tfid = []
i = -1
gene =[]

'''
读出具有一一对应关系的DBD与family名称
'''
family_id = []
DBD = []
with open(family_db,'r') as db:
    for line in db:
        line = line.strip()
        line = line.split("\t")
        family_id.append(line[0])
        DBD.append(line[1])


'''
逐行分析比对结果
提取输入序列的名称
& 对应的比对内容
'''
for line in out:
    line = line.strip()
    tem = re.match('#.*', line)
    if tem:
        pass
    elif line == '':
        pass
    elif 'Query:' in line:
        x = line.split(' ')
        query.append(x[7])
        i += 1
        phmmer.append("")
    elif 'Description:' in line:
        line = line.strip()
        g = re.search('gene:([a-zA-Z0-9\.]*)\s', line)
        gene.append(g.group(1))
        i += 1
        phmmer.append("")
    else:
        phmmer[i] += line+"\n"

'''
按照文献中的方式区分家族
'''
family = []
for i in range(len(query)):
    x = re.findall(">>\s(.*)\s{2}", phmmer[i])
    if not x:
        family.append("None")
    elif x[0].lstrip() in DBD:
        inx = DBD.index(x[0].lstrip())
        family.append(family_id[inx]) 
    elif "AP2" in x:
        if "B3" in x:
            family.append("RAV")
        elif len(x) >= 2:
            family.append("AP2")
        else:
            family.append("ERF")
    elif "B3" in x:
        if "Auxin_resp" in x:
            family.append("ARF")
        else:
            family.append("B3")
    elif "zf-B_box" in x:
        if "CCT" in x:
            family.append("CO-like")
        elif len(x) >= 2:
            family.append("DBB")
        else:
            family.append("None")
    elif "Zf-LSD1" in x:
        if "Peptidase_C14" in x:
            family.append("None")
        else:
            family.append("LSD")
    elif "zf-C2H2" in x:
        if "RNase_T" in x:
            family.append("None")
        else:
            family.append("C2H2")
    elif "Zf-CCCH" in x:
        if "RRM_1" in x or "Helicase_C" in x:
            family.append("None")
        else:
            family.append("C3H")
            family.append("C3H")
    elif "G2-like" in x:
        if "Response_reg" in x:
            family.append("ARR-B")
        else:
            family.append("G2-like")
    elif "WRC" in x and "QLQ" in x:
        family.append("GRF")
    elif "Homeobox" in x:
        if "HD-ZIP_I/II" in x or "SMART" in x:
            family.append("HD-ZIP")
        elif "BELL" in x or "ELK" in x:
            family.append("TALE")
        elif "Wus type homeobox" in x:
            family.append("WOX")
        elif "PHD" in x:
            family.append("HB-PHD")
        else:
            family.append("HB-other")
    elif "SRF-TF" in x:
        if "K-box" in x:
            family.append("MIKC")
        else:
            family.append("M_type")
    elif "Myb_dna_bind" in x:
        if len(x) >= 2 and "SWIRM" not in x:
            family.append("MYB")
        elif len(x) == 1 and "SWIRM" not in x:
            family.append("MYB_related")
    else:
        family.append("None")

diff_gene_id = sys.argv[4]
# diff_gene_id = 'featurecounts_count.txt.A_vs_B.edgeR.DE_results_name'
diff = []
with open(diff_gene_id, 'r') as d:
    for line in d:
        line = line.strip()
        diff.append(line)
print len(family)
print len(query)

for k in range(len(family)):
    if family[k] == "None":
        pass
    else:
        if gene[k] in diff:
            # print query[k] +'\t'+ family[k] +'\t'+gene[k] +'\n'
            f_out = open("TF_result.txt","a+")
            f_out.write(query[k] +'\t'+ family[k] +'\t'+gene[k] +'\n')
            f_out.close()
        else:
            # print query[k] +'\t'+ family[k] +'\t'+'no' +'\n'
            f_out = open("TF_result.txt","a+")
            f_out.write(query[k] +'\t'+ family[k] +'\t'+'no' +'\n')
            f_out.close()
