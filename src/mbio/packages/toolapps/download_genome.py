#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# 根据GCF和GCA号下载对应的基因数据
import os
import re
import argparse
import gzip
import time
from Bio import SeqIO

parser = argparse.ArgumentParser(description="download genome data from NCBI")
parser.add_argument('-ID_list', required=True, type=str, help='requested genome number list')
parser.add_argument('-output', required=True, type=str, help='output_dir')
parser.add_argument('-path_list', required=False, default="file.txt", type=str, help='path list of gff/fna/report')
parser.add_argument('-sh', required=False, default="download.sh", type=str, help='download sh file')
args = parser.parse_args()
args.output = os.path.abspath(str(args.output))
args.ID_list = str(args.ID_list)
args.path_list = str(args.path_list)
args.sh = str(args.sh)
if not os.path.exists(args.output):
    os.mkdir(args.output)

###根据GBK文件提取物种分类信息，并将分类信息和物种信息连起来
def get_taxon(file):
   species = ""
   ta = ""
   for rec in SeqIO.parse(file, "gb"):
      for key in rec.annotations.keys():
          if key in ['organism']:
             species = rec.annotations['organism']
          if key	in ['taxonomy']:
             ta = ";".join(rec.annotations['taxonomy'])
   return ta+";"+species

###根据rna文件获取16s的个数以及最长16s的序列文件，并改名字
def get_rna_fasta(fa,sample,out):
   num =0
   list=[]
   for i in SeqIO.parse(fa, "fasta"):
      if re.search("product=16S ribosomal RNA",i.description):
        num +=1
        i.id = sample+"_rRNA"+str(num)
        list.append(i)
   SeqIO.write(list, out, "fasta")
   return num

fi = open(args.sh, "w")
gca_gcf = args.output + "/gca_gcf.txt"
f_gcf = open(gca_gcf, "a")
list = args.ID_list.split(";")
for name in list:
   pr2 = ''
   if 'GCF' in name:
       res1 = os.popen("esearch -db assembly -query " + name + " | efetch -format docsum ").read()
       if res1:
          pr = re.findall("<FtpPath_RefSeq>(.*?)<", res1)
          pr1 = pr[0].split('/', 3)[3]
          pr2 = "ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k1 -QTr -l100m anonftp@ftp.ncbi.nlm.nih.gov:" + pr1 + " " + args.output
   elif 'GCA' in name:
       res1 = os.popen("esearch -db assembly -query " + name + " | efetch -format docsum ").read()
       if res1:
          pr = re.findall("<FtpPath_GenBank>(.*?)<", res1)
          if pr:
              pr1 = pr[0].split('/', 3)[3]
          else:
              pr_tmp = re.findall("<FtpPath_RefSeq>(.*?)<", res1)
              pr1 = pr_tmp[0].split('/', 3)[3]
              gcf = re.findall("<RefSeq>(.*?)<", res1)
              f_gcf.write("{}\t{}\n".format(name, "".join(gcf)))
          pr2 = "ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k1 -QTr -l100m anonftp@ftp.ncbi.nlm.nih.gov:" + pr1 + " " + args.output

   else:
       print(name + " is not true genome number")
   if pr2:
       fi.write(pr2 + "\n")
   time.sleep(30)
fi.close()
f_gcf.close()
os.system("sh" + " " + args.sh)

# 生成path_list文件，主要存放gff，fna以及report的绝对路径信息，为了后边的处理, 还有解压gz文件的功能
path = open(args.path_list, "w")
for i in os.listdir(args.output):
    if i == "gca_gcf.txt":
        pass
    else:
        d1 = os.path.join(args.output, i)
        full_name = d1.split("/")[-1]
        sim_name = full_name.split("_")[0] + "_" + full_name.split("_")[1]
        fna_path = d1 + "/" + full_name + "_genomic.fna.gz"
        fna = fna_path.replace(".gz", "")
        f = gzip.open(fna_path, "rb")
        f_content = f.read()
        f.close()
        open(fna, "wb+").write(f_content)
        fna_path_1 = d1 + "/" + full_name + "_genomic.fna"
        gbk_path = d1 + "/" + full_name + "_genomic.gbff.gz"
        taxon = ''
        if not os.path.exists(gbk_path):
            gbk_path_1 = "-"
        else:
            gbk = gbk_path.replace(".gz", "")
            g = gzip.open(gbk_path, "rb")
            g_content = g.read()
            g.close()
            open(gbk, "wb+").write(g_content)
            gbk_path_1 = d1 + "/" + full_name + "_genomic.gbff"
            taxon = get_taxon(gbk_path_1)
        rna_path = d1 + "/" + full_name + "_rna_from_genomic.fna.gz"
        num = 0
        rna_path_2 = d1 + "/" + full_name + "_rrna.fna"
        if not os.path.exists(rna_path):
            rna_path_2 = "-"
        else:
            rna = rna_path.replace(".gz", "")
            f = gzip.open(rna_path, "rb")
            f_content = f.read()
            f.close()
            open(rna, "wb+").write(f_content)
            rna_path_1 = d1 + "/" + full_name + "_rna_from_genomic.fna"
            rna_path_2 = d1 + "/" + full_name + "_rrna.fna"
            num = get_rna_fasta(rna_path_1, sim_name, rna_path_2)
        input_line = sim_name + "\t" + fna_path_1 + "\t" + taxon + "\t" + str(num) + "\t" + rna_path_2+"\n"
        path.write(input_line)
path.close()


