# -*- coding: utf-8 -*-
import os
import argparse
from collections import defaultdict
parser = argparse.ArgumentParser()
parser.add_argument("-new",type=str,required=True,help="please input correct before stat info")
parser.add_argument("-known",type=str,required=True,help="please input correct after stat info")
parser.add_argument("-newm",type=str,required=True,help="please input correct before stat info")
parser.add_argument("-knownm",type=str,required=True,help="please input correct after stat info")
parser.add_argument("-quant",type=str,required=True,help="please input quant file")
parser.add_argument("-out",type=str,required=True,help="please input correct out name")
parser.add_argument("-t2g",type=str,help="please input correct out name")#非必需，输入文件为gene相关的的时候需要
args=parser.parse_args()

new_lnc=args.new
known_lnc=args.known
new_m=args.newm
known_m=args.knownm
out_dir=args.out
quant_dir=args.quant
lnc_rna_list=[]
mrna_list=[]
new_rna_list=[]
ref_rna_list=[]


with open(new_lnc,"r")as new_lnc_ids:
    for new_lnc_id in new_lnc_ids.readlines():
        new_lnc_id=new_lnc_id.strip().split()[0]
        lnc_rna_list.append(new_lnc_id)
        new_rna_list.append(new_lnc_id)
with open(known_lnc, "r")as known_lnc_ids:
    for known_lnc_id in known_lnc_ids.readlines():
        known_lnc_id=known_lnc_id.strip().split()[0]
        lnc_rna_list.append(known_lnc_id)
        ref_rna_list.append(known_lnc_id)
with open(new_m, "r")as new_m_ids:
    for new_m_id in new_m_ids.readlines():
        new_m_id=new_m_id.strip().split()[0]
        mrna_list.append(new_m_id)
        new_rna_list.append(new_m_id)
with open(known_m, "r")as known_m_ids:
    for known_m_id in known_m_ids.readlines():
        known_m_id=known_m_id.strip().split()[0]
        mrna_list.append(known_m_id)
        ref_rna_list.append(known_m_id)


if args.t2g :
    t2g=args.t2g
    g2t=defaultdict(list)
    with open(args.t2g,"r") as t2g_infos:
        for t2g_info in t2g_infos.readlines():
            g2t[t2g_info.strip().split()[1]].append(t2g_info.strip().split()[0])
    with open(out_dir,"w")as final_quant_results:
      with open(quant_dir,"r")as quant_resluts:
        first_info_quant=quant_resluts.readline().strip()
        final_quant_results.write(first_info_quant+"\t"+"rna_type"+"\t"+"is_new"+"\n")
        for quant_detail_resluts in quant_resluts.readlines():
            gene_id=quant_detail_resluts.strip().split()[0]
            trans_id=g2t[gene_id]
            l=0
            m=0
            new=0
            ref=0
            for tr_id in trans_id:
                if tr_id in lnc_rna_list:
                    l+=1
                elif tr_id in mrna_list:
                    m+=1
                else:
                    continue
            for tr_id in trans_id:
                if tr_id in new_rna_list:
                    new+=1
                elif tr_id in ref_rna_list:
                    ref+=1
                else:
                    continue
            if m>=1:
                if ref>=1:
                    final_quant_results.write(quant_detail_resluts.strip() + "\t" + "mRNA" + "\t"+"false"+"\n")
                elif new>=1 and ref==0:
                    final_quant_results.write(quant_detail_resluts.strip() + "\t" + "mRNA" + "\t" + "true" + "\n")
            elif l>=1 and m==0:
                if ref>=1:
                     final_quant_results.write(quant_detail_resluts.strip() + "\t" + "lncRNA"  + "\t"+"false"+"\n")
                elif new>=1 and ref==0:
                    final_quant_results.write(quant_detail_resluts.strip() + "\t" + "lncRNA" + "\t" + "true" + "\n")
            else:
                continue
    print(g2t)

else:
    with open(out_dir, "w")as final_quant_results:
        with open(quant_dir, "r")as quant_resluts:
            first_info_quant = quant_resluts.readline().strip()
            final_quant_results.write(first_info_quant + "\t" + "rna_type" +"\t"+"is_new"+"\n")
            for quant_detail_resluts in quant_resluts.readlines():
                trans_id = quant_detail_resluts.strip().split()[0]
                if trans_id in lnc_rna_list:
                    if trans_id in ref_rna_list:
                            final_quant_results.write(quant_detail_resluts.strip() + "\t" + "lncRNA" + "\t"+"false"+"\n")
                    else:
                            final_quant_results.write(quant_detail_resluts.strip() + "\t" + "lncRNA" + "\t" + "true" + "\n")
                elif trans_id in mrna_list:
                    if trans_id in ref_rna_list:
                             final_quant_results.write(quant_detail_resluts.strip() + "\t" + "mRNA" +"\t"+"false"+"\n")
                    else:
                             final_quant_results.write(quant_detail_resluts.strip() + "\t" + "mRNA" + "\t" + "true" + "\n")
                else:
                    continue



