# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171208

"""对比对到nt数据库的table文件进行统计，得到物种的统计表"""
import argparse


parser = argparse.ArgumentParser(description='Compare the table results to the nt database for species statistics')
parser.add_argument('-i', '--table', help='input table file', required=True)
parser.add_argument('-f', '--fasta', help='input fasta file', required=True)
parser.add_argument('-n', '--number', help='input species number,  default:5', required=False)
parser.add_argument('-o', '--output', help='output stat file', required=True)
args = vars(parser.parse_args())

if args['number']:
    number = int(args['number'])
else:
    number = 5
table_file = args['table']
fasta_file = args['fasta']
outfile = args['output']

sp = {}
with open(table_file, "r") as f:
    lines = f.readlines()
    blast_sum = len(lines) - 1
    for line in lines[1:]:
        item = line.strip().split("\t")
        des = item[15].split()
        if des[1] == "PREDICTED:":
            # species = des[1:3]
            species = des[2:4]
        else:
            # species = des[0:2]
            species = des[1:3]
        species = " ".join(species)
        if species not in sp:
            sp[species] = 1
        else:
            sp[species] += 1
fa = open(fasta_file)
fa_sum = len(fa.readlines()) / 2
with open(outfile, "w") as w:
    w.write("Species\tNum\tPercent(%)\n")
    sp_ = sorted(sp.iteritems(), key=lambda item: item[1], reverse=True)
    other = 0
    for i in range(len(sp_)):
        if i <= number - 1:
            line = []
            line.append(sp_[i][0])
            line.append(str(sp_[i][1]))
            line.append(str(round(float(sp_[i][1]) / fa_sum * 100, 2)))
            print line
            w.write("\t".join(line) + "\n")
        else:
            other += sp_[i][1]
    if other != 0:
        line = []
        line.append("Other species")
        line.append(str(other))
        line.append(str(round(float(other) / fa_sum * 100, 2)))
        w.write("\t".join(line) + "\n")
    line = []
    noblast = fa_sum - blast_sum
    line.append("No blast")
    line.append(str(noblast))
    line.append(str(round(float(noblast) / fa_sum * 100, 2)))
    w.write("\t".join(line) + "\n")
