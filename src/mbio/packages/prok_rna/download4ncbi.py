#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# 根据GCF和GCA号下载对应的基因数据
import os
import re
import argparse
import gzip
import time
from BCBio import GFF

parser = argparse.ArgumentParser(description="download genome data from NCBI")
parser.add_argument('-ID_list', required=True, type=str, help='requested genome number list')
parser.add_argument('-output', required=True, type=str, help='output_dir')
parser.add_argument('-path_list', required=False, default="file.txt", type=str, help='path list of gff/fna/report')
parser.add_argument('-sh', required=False, default="download.sh", type=str, help='download sh file')
parser.add_argument('-sorfware', required=True, type=str, help='esearch path')
args = parser.parse_args()
args.output = os.path.abspath(str(args.output))
args.ID_list = str(args.ID_list)
args.path_list = str(args.path_list)
args.sh = str(args.sh)
if not os.path.exists(args.output):
    os.mkdir(args.output)

fi = open(args.sh, "w")
with open(args.ID_list, "r") as list:
    for line in list:
        name = line.strip("\n")
        if 'GCF' in name:
            res1 = os.popen("{}/bioinfo/ref_rna_v3/gene_fusion/miniconda3/bin/esearch -db assembly -query ".format(args.sorfware) + name + " | efetch -format docsum ").read()
            pr = re.findall("<FtpPath_RefSeq>(.*?)<", res1)
            pr1 = pr[0].split('/', 3)[3]
            pr2 = "~/.aspera/connect/bin/ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k1 -QTr -l100m anonftp@ftp.ncbi.nlm.nih.gov:" + pr1 + " " + args.output
        elif 'GCA' in name:
            res1 = os.popen("{}/bioinfo/ref_rna_v3/gene_fusion/miniconda3/bin/esearch -db assembly -query ".format(args.sorfware) + name + " | efetch -format docsum ").read()
            pr = re.findall("<FtpPath_GenBank>(.*?)<", res1)
            pr1 = pr[0].split('/', 3)[3]
            pr2 = "~/.aspera/connect/bin/ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k1 -QTr -l100m anonftp@ftp.ncbi.nlm.nih.gov:" + pr1 + " " + args.output
        else:
            print(name + " is not true genome number")
        fi.write(pr2 + "\n")
        time.sleep(30)
    fi.close()
    os.system("sh" + " " + args.sh)

# 生成path_list文件，主要存放gff，fna以及report的绝对路径信息，为了后边的处理, 还有解压gz文件的功能
path = open(args.path_list, "w")
for i in os.listdir(args.output):
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
    gff_path = d1 + "/" + full_name + "_genomic.gff.gz"
    gtf_path = d1 + "/" + full_name + "_genomic.gtf.gz"
    if not os.path.exists(gff_path):
        if os.path.exists(gtf_path):
            gtf = gtf_path.replace(".gz", "")
            g = gzip.open(gtf_path, "rb")
            g_content = g.read()
            g.close()
            open(gtf, "wb+").write(g_content)
            gff_path_1 = d1 + "/" + full_name + "_genomic.gtf"
        else:
            gff_path_1 = '-'
    else:
        gff = gff_path.replace(".gz", "")
        g = gzip.open(gff_path, "rb")
        g_content = g.read()
        g.close()
        open(gff, "wb+").write(g_content)
        gff_path_1 = d1 + "/" + full_name + "_genomic.gff"
        # 判断gff文件格式是否通过检查
        in_handle = open(gff_path_1)
        try:
            for rec in GFF.parse(in_handle):
                pass
        except:
            if os.path.exists(gtf_path):
                gtf = gtf_path.replace(".gz", "")
                g = gzip.open(gtf_path, "rb")
                g_content = g.read()
                g.close()
                open(gtf, "wb+").write(g_content)
                gff_path_1 = d1 + "/" + full_name + "_genomic.gtf"
            else:
                gff_path_1 = '-'
    report_path = d1 + "/" + full_name + "_assembly_report.txt"
    input_line = sim_name + "\t" + fna_path_1 + "\t" + gff_path_1 + "\t" + report_path
    path.write(input_line + "\n")
path.close()
