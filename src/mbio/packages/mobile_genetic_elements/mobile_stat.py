# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from optparse import OptionParser
from Bio import SeqIO
import re,os
import pandas as pd

def change_result(genome, dir, name, prefix):
    mge = {}
    elem = {}
    files = os.listdir(dir)
    m=1
    for file in files:
        if os.path.getsize(dir + "/" + file) >0:
            if re.search('integrons', file):
                with open(dir + "/" + file, "r") as f:
                    lines = f.readlines()
                    integron = {}
                    for line in lines:
                        if re.search("^#",line):
                            continue
                        elif re.search("^ID_integron",line):
                            continue
                        else:
                            lin = line.strip().split("\t")
                            des = lin[0] + "\t" + lin[1]
                            if des not in integron.keys():
                                integron[des] = [lin]
                            else:
                                integron[des].append(lin)
                    for key in integron.keys():
                        postion = []
                        location = ''
                        names = ''
                        for i in integron[key]:
                            location = i[1]
                            names = i[13]
                            postion.append(int(i[3]))
                            postion.append(int(i[4]))
                            if i[8] == "attC":
                                length1 = abs(int(i[3]) - int(i[4])) + 1
                                if i[5] == '1':
                                    stran = "+"
                                    elem["elem" + str(m)] = [i[8], i[9], i[1], int(i[3]), i[4], stran, length1]
                                elif i[5] == '-1':
                                    stran = "-"
                                    elem["elem" + str(m)] = [i[8], i[9], i[1], int(i[4]), i[3], stran, length1]
                                m +=1
                        postion = sorted(postion)
                        length = abs(int(postion[-1])-int(postion[0])) +1
                        mge[key] =["integron", "integron", names, location, int(postion[0]), postion[-1], '+', length, "-"]
            elif re.search('island', file):
                with open(dir + "/" + file, "r") as f:
                    lines = f.readlines()
                    for line in lines:
                        if re.search("^Location",line):
                            continue
                        else:
                            lin = line.strip().split("\t")
                            mge[lin[1]] = ["island", "island","island",lin[0], int(lin[2]), lin[3], '+', lin[4], lin[5]]
            elif re.search('transposon', file):
                with open(dir + "/" + file, "r") as f:
                    lines = f.readlines()
                    n = 1
                    for line in lines[1:]:
                        lin = line.strip().split("\t")
                        if lin[0] == "transposon":
                            mge["transposon"+str(n)] =["transposon", "transposon","transposon", lin[2], int(lin[3]), lin[4], '+', lin[5], "-"]
                            n += 1
                        elif lin[0] == "is":
                            length3 = abs(int(lin[3]) - int(lin[4])) + 1
                            mge["transposon" + str(n)] = ["is", lin[6], lin[7], lin[2], int(lin[3]), lin[4], lin[5], length3, "-"]
                            if lin[8] == "-":
                                pass
                            else:
                                ir1 = lin[8].split("..")
                                ir2 = lin[9].split("..")
                                length1 = abs(int(ir1[0]) - int(ir1[1])) + 1
                                elem["elem" + str(m)] = ["IR", "IR1", lin[2], int(ir1[0]), ir1[1], "+", length1]
                                m += 1
                                length2 = abs(int(ir2[0]) - int(ir2[1])) + 1
                                elem["elem" + str(m)] = ["IR", "IR2", lin[2], int(ir2[0]), ir2[1], "+", length2]
                                m += 1
                        n += 1
            elif re.search('prompter', file):
                with open(dir + "/" + file, "r") as f:
                    lines = f.readlines()
                    n = 1
                    for line in lines[1:]:
                        lin = line.strip().split("\t")
                        mge["prompter"+str(n)] = ['prompter', 'prompter', 'prompter', lin[0], int(lin[1]), lin[2], lin[3], lin[4], "-"]
                        n += 1
            elif re.search('enzyme', file):
                with open(dir + "/" + file, "r") as f:
                    lines = f.readlines()
                    for line in lines:
                        lin = line.strip().split("\t")
                        length = abs(int(lin[2]) - int(lin[3]))+ 1
                        elem["elem"+str(m)] = ['enzyme', lin[5], lin[1], int(lin[2]), lin[3], lin[4], length]
                        m += 1
    ##处理MGE的数据
    if len(mge.keys()) >0:
        dict_mge = get_seq(genome, mge, 3, 4, 5, 6)
        for key in mge.keys():
            if key in dict_mge.keys():
                mge[key].append(dict_mge[key])
        data = pd.DataFrame.from_dict(mge).T
        data.columns = ['type', "mge_group", "mge_name", "location", "start", "end", "strand", "length", "sofware",
                        "seq"]
        data['start'].astype('int')
        data['data_type'] = "mge"
        data['sample'] = name
        data.sort_values(["location", "start"], inplace=True)
        data.to_csv(prefix + "/" + name + ".mge.xls", sep='\t', header=True, index=False,
                    columns=['sample', 'data_type', 'type', "mge_group", "mge_name", "location", "start", "end",
                             "strand", "length", "sofware", "seq"])
    ##处理element的数据
    if len(elem.keys()) >0:
        dict_elem = get_seq(genome, elem, 2, 3, 4, 5)
        for key in elem.keys():
            if key in dict_elem.keys():
                elem[key].append(dict_elem[key])
        data2 = pd.DataFrame.from_dict(elem).T
        data2.columns = ['type', "elem_name", "location", "start", "end", "strand", "length", "seq"]
        data2['data_type'] = "elem"
        data2['sample'] = name
        data2.sort_values(["location", "start"], inplace=True)
        data2.to_csv(prefix + "/" + name + ".element.xls", sep='\t', header=True, index=False,
                     columns=['sample', 'data_type', 'type', "elem_name", "location", "start", "end", "strand",
                              "length", "seq"])

def get_seq(fa, dict1, loction, start, end, strand):
    dict_seq = {}
    for key in dict1.keys():
        for seq_record in SeqIO.parse(fa, "fasta"):
            if dict1[key][loction] == seq_record.id:
                if dict1[key][strand] == "+":
                    seq = str(seq_record.seq[int(dict1[key][start]) - 1:int(dict1[key][end])])
                    dict_seq[key] = seq
                elif dict1[key][strand] == "-":
                    seq1 = str(seq_record.seq[int(dict1[key][end]) - 1:int(dict1[key][start])])
                    seq = dna_reverse(dna_complement(seq1))
                    dict_seq[key] = seq
    return dict_seq

def dna_complement(sequence):
    sequence = sequence.upper()
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')
    return sequence.upper()

def dna_reverse(sequence):
    sequence = sequence.upper()
    return sequence[::-1]

def main():
    parser = OptionParser()
    parser.add_option('--n', dest='name', metavar='[sample name]')
    parser.add_option('--g', dest='genome', metavar='[genome fna file]')
    parser.add_option('--d', dest='dir', metavar='[files dir]')
    parser.add_option("--o", dest="prefix", metavar="[prefix of  file]")
    (options,args) = parser.parse_args()
    if not options.name or not options.genome or not options.dir or not options.prefix:
        print "python mobile_stat.py --n name --d dir --o prefix"
        return
    change_result(options.genome, options.dir, options.name, options.prefix)

if __name__=='__main__':
    main()