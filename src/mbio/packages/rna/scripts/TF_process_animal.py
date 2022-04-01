# -*- coding: utf-8 -*-
# author: Moli
# last update: 20161108

#python TF_process_animal_v1.py ../hmmerdb/animal/animaltfdb.hmm
# /mnt/ilustre/users/sanger-dev/app/database/refGenome/Animal/Primates/human/Homo_sapiens.GRCh38.pep.all.fa
# family_vs_DBD_animal_2.0.txt
# /mnt/ilustre/users/sanger-dev/workspace/20161116/Single_express_sample_v8_change_sample_name
# /Express/output/diff/featurecounts_count.txt.A_vs_B.edgeR.DE_results_name

#python TF_process_animal_v2.py family_vs_DBD_animal_2.0.txt /mnt/ilustre/users/sanger-dev/workspace/20161116/Single_express_sample_v8_change_sample_name/Express/output/diff/featurecounts_count.txt.A_vs_B.edgeR.DE_results_name

import re
import sys
import subprocess
#
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
#
family_db = sys.argv[3]
# family_db = 'O:/Users/moli.zhou/Desktop/family_vs_DBD_animal_2.0.txt'
query = [] #输入的蛋白序列号
phmmer = [] #hmmscan比对结果
domain = []
tfid = []
gene = [] #蛋白对应的gene序列号
i = -1

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
# file = 'Homo_sapiens_hmm.out'
# with open(file, 'r') as f:
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


family =[]
for i in range(len(phmmer)):
    if ">>" in phmmer[i]:
        y = re.search(">>\s(.*)\s", phmmer[i])
        z = y.group(1).split()
        for j in range(len(DBD)):
            if z[0] == DBD[j]:
                family.append(family_id[j])
        if z[0] not in DBD:
            family.append("Other")
    else:
        family.append("None")

diff_gene_id = sys.argv[4]
# diff_gene_id = 'featurecounts_count.txt.A_vs_B.edgeR.DE_results_name'
diff = []
with open(diff_gene_id, 'r') as d:
    for line in d:
        line = line.strip()
        diff.append(line)

for k in range(len(query)):
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
