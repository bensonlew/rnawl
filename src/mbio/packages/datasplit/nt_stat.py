# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171208

"""对比对到nt数据库的table文件进行统计，得到物种的统计表"""
import argparse
from biocluster.config import Config


parser = argparse.ArgumentParser(description='Compare the table results to the nt database for species statistics')
parser.add_argument('-i', '--table', help='input table file', required=True)
parser.add_argument('-f', '--fasta', help='input fasta file', required=True)
parser.add_argument('-n', '--number', help='input species number,  default:5', required=False)
parser.add_argument('-o', '--output', help='output stat file', required=True)
args = vars(parser.parse_args())

# coll_tax = Config().biodb_mongo_client.sanger_biodb.NT_sequence_20200604
# coll_spe = Config().biodb_mongo_client.sanger_biodb.NT_taxid2species
coll_tax = Config().biodb_mongo_client.sanger_biodb1.NT_sequence_20200604
coll_spe = Config().biodb_mongo_client.sanger_biodb1.NT_taxid2species20211009
# coll_spe = Config().biodb_mongo_client.sanger_biodb1.NT_taxid2species
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
        # acci_id = item[15].split()[0]
        acci_id = item[4].split("|")[-2]
        result = coll_tax.find_one({"_id": acci_id})
        if result:
            taxid = result["taxid"]
            result1 = coll_spe.find_one({"_id": int(taxid)})
            if result1:
                species = result1["species"]
                if species not in sp:
                    sp[species] = 1
                else:
                    sp[species] += 1
            else:
                print "NT_taxid2species里没有找到_id:"+taxid
        else:
            print "NT_sequence_20200604里没有找到_id:"+acci_id
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
