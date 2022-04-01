#/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = '刘彬旭'

"""
genome json
"""

import sys
import json
from optparse import OptionParser
reload(sys)
sys.setdefaultencoding('utf8')
import os

table_file = sys.argv[1]
json_file = table_file.split(".xls")[0] + ".json"

genome_dict = dict()
with open(table_file, 'r') as f:
    head = f.readline().strip("\n").strip("#").split("\t")
    line = f.readline().strip("\n").strip("#").split("\t")
    genome_dict = dict(zip(head, line))

print genome_dict

if genome_dict['ensemble_class'] not in ['protists', 'metazoa', 'vertebrates', 'plants', 'fungi']:
    raise Exception("{} 不能为空且必须为以下之一\n{}".format('ensemble_class', ['protists', 'metazoa', 'vertebrates', 'plants', 'fungi']))
if genome_dict['taxon'] not in ['Animal', 'Fungi', 'Plant', 'Protist', 'fungi']:
    raise Exception("{} 不能为空且必须为以下之一\n{}".format('taxon', ['Animal', 'Fungi', 'Plant', 'Protist', 'fungi']))
if genome_dict['ensemble_release'] == '':
    raise Exception("ensemble_release 不能为空, 物种数据所在目录")

if genome_dict['name'] == '':
    raise Exception("name 不能为空, 页面显示的物种名")
if genome_dict['index'] == '':
    raise Exception("index 不能为空, 页面显示的")
if genome_dict['biomart_gene_annotype'] == '':
    print "注意 biomart_gene_annotype 为空， 默认设为type3"
if genome_dict['assembly'] == '':
    print "注意 assembly 为空， 默认用ensemble_release + accession两个字段决定, 该字段会在参数页面选择"
if genome_dict['annot_version'] == '':
    print "注意 annot_version 为空， 默认用ensemble_release 字段决定， 该字段会在参数页面选择"



abs_path = os.path.abspath(table_file)
now_dir = os.getcwd()
dirs = abs_path.split("/")

if dirs[-2] !=  genome_dict['ensemble_release']:
    print '不在 {} 目录下, 将文件转移进去或者检查表格中的ensemble_release 字段'.format(genome_dict['ensemble_release'])
if dirs[-3] !=  genome_dict['name']:
    print '不在 {} 目录下, 将文件转移进去或者检查表格中的name字段'.format(genome_dict['name'])
if dirs[-4] !=  genome_dict['ensemble_class']:
    print '不在 {} 目录下, 将文件转移进去或者检查表格中的ensemble_class字段'.format(genome_dict['ensemble_class'])

genome_base_dir = "/".join(dirs[:-4])
genome_dir = "/".join(dirs[-4:-1])

files = [
    'ensemble2entrez',
    'dna_fa',
    'dna_index',
    'cds',
    'pep',
    'bio_mart_annot',
    'go',
    'kegg',
    'cog',
    'gtf',
    'gene_stat',
    # 'anno_path',
    # 'anno_path_v2',
]
afile2abr = {
    'ensemble2entrez': ".enterz.txt",
    'dna_fa': ".fa",
    'dna_index': "",
    'cds': ".cds.fa",
    'pep': ".pep.fa",
    'bio_mart_annot': ".biomart",
    'go': ".gene2go",
    'kegg': ".pathway",
    'cog': ".gene2cog.xls",
    'gtf': ".gtf",
    'gene_stat': ".genome_stat.xls",
    # 'anno_path',
    # 'anno_path_v2':,
}

afile2dir = {
    'ensemble2entrez': "NCBI",
    'dna_fa': "dna",
    'dna_index': "dna",
    'cds': "cds",
    'pep': "cds",
    'bio_mart_annot': "biomart",
    'go': "GO",
    'kegg': "KEGG",
    'cog': "COG",
    'gtf': "gtf",
    'gene_stat': "gtf",
    # 'anno_path',
    # 'anno_path_v2':,
}
# check file
def check_file(afile, file_path):
    # print file_path, genome_dir
    if not file_path.startswith(genome_dir):
        print 'warning {} file not in {} are you sure?'.format(file_path, genome_dir)
    file_path = genome_base_dir + '/' + file_path
    if afile == 'dna_index':
        file_path += '.1.ht2'
        file_path2 = file_path + "l"
        if os.path.exists(file_path) or os.path.exists(file_path2):
            pass
        else:
            print 'warning {} file not exists'.format(genome_dict[afile])
    else:
        if os.path.exists(file_path):
            pass
        else:
            print 'warning {} file not exists'.format(genome_dict[afile])

for afile in files:
    if genome_dict[afile] != "":
        file_path = genome_dict[afile]
        check_file(afile, file_path)
    else:
        genome_dict[afile]  = genome_dir + "/" + afile2dir[afile] + "/" + genome_dict['index'] + afile2abr[afile]
        file_path = genome_dict[afile]
        check_file(afile, file_path)

with open(genome_base_dir + '/' + genome_dict['gene_stat'], 'r') as f:
    f.readline()
    size = 0
    for line in f.readlines():
        size += float(line.split("\t")[1])

if genome_dict['size'] == "":
    genome_dict['size'] = '%.2f' % size

if genome_dict['biomart_gene_annotype'] == "":
    genome_dict['biomart_gene_annotype'] = "type3"
if genome_dict['taxon_id'] == "":
    genome_dict['taxon_id'] = 0

genome_dict["organism_name"] = genome_dict["name"]
genome_dict["anno_path_v2"] = genome_dir + "/Annotation_v2"

# 推测基因组版本号
if genome_dict.has_key('genome_version') and genome_dict['genome_version']:
    genome_dict['assembly'] = genome_dict['genome_version']
elif genome_dict['assembly'] != "":
    pass
elif genome_dict['ensemble_release'] != "" or genome_dict['accession'] != "":
    genome_dict['assembly'] = genome_dict['ensemble_release'] + "_" + genome_dict['accession']
else:
    genome_dict['assembly'] = "unknown"

# 推测注释版本号
if genome_dict['annot_version']:
    genome_dict['annot_version'] = genome_dict['annot_version']
elif genome_dict['ensemble_release'] != "":
    genome_dict['annot_version'] = genome_dict['ensemble_release']
else:
    genome_dict['annot_version'] = "i-sanger"
# print genome_dict


jsonout = open(json_file,'w')
jsonout.write(json.dumps(genome_dict, indent=4))
jsonout.close()

