# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

import argparse
import os
import subprocess
from turtle import tilt, tiltangle

parser = argparse.ArgumentParser(description='''''')
parser.add_argument("-i","--input", help="输入blast后的结果文件，在-outfmt \"6 -delim \"\t\" qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids \"", required=True)
parser.add_argument("-o", "--out", help="输入输出路径", required=True)
args = vars(parser.parse_args())


In_NT = args["input"]
out_NT = args["out"]
taxonkit = "/mnt/ilustre/users/sanger-dev/app/bioinfo/tool_lab/taxonkit"
# temp = "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/NT/test/temp"

def getTaxon(taxid,is_print):
    code = subprocess.Popen("echo {}|{} lineage|{} reformat -f \"{{k}};{{p}};{{c}};{{o}};{{f}};{{g}};{{s}}\" -P -r norank|cut -f 3".format(taxid,taxonkit,taxonkit),shell = True, stdout = subprocess.PIPE,
            stderr = subprocess.STDOUT)
    taxname = ""
    code.wait()
    if code.poll() != 0:
        raise Exception("获取Taxon程序运行失败，返回code:{}".format(code.stdout.read()))
    # print code.stdout.read()
    try:
        taxname = code.stdout.read()
        # print taxname
        # with open(temp,"r") as t:
        #     taxname = t.readline()
        if is_print:
            print taxname
    except Exception as e:
        print("taxon名没有获取到{}".format(taxid))
    # if taxname.strip() == "":
    #     print taxid
    return taxname.strip()

with open(In_NT,"r") as inref, open(out_NT,"w") as out_nt:
    old_taxid = ""
    is_show = False
    old_taxname = ""
    while 1:
        line = inref.readline()
        if not line:
            break
        fd = line.rstrip().split('\t')
        taxid = fd[-1]
        accession_id = fd[1].split('|')[3]
        if taxid == old_taxid:
            taxname = old_taxname
        else:
            if len(taxid.split(';')) >1:
                taxname_last = ""
                taxnames = []
                # is_show = True
                # print taxid
                for i in taxid.split(';'):
                    taxname_tem = getTaxon(i,is_show)
                    if taxname_last != taxname_tem:
                        taxname_last = taxname_tem
                        taxnames.append(taxname_tem)
                taxname = ";".join(taxnames)
                old_taxid = taxid
                old_taxname = taxname
            else:
                taxname = getTaxon(taxid,is_show)
                is_show = False
                old_taxid = taxid
                old_taxname = taxname
        out_nt.write("{}\t{}\t{}\t".format(fd[0],taxname.replace(' ','_'),accession_id))
        out_nt.write("\t".join(fd[2:-1])+"\n")



            


