#!/usr/bin/env python
# -*- coding: utf-8 -*-
# _auther_ zhenyu.lu
# _time_ 2020/6/22


import os
import argparse
import subprocess
import sys
import re
import pandas as pd
import numpy as np


def arrange_gtf(gtfile,outdir):
    with open(gtfile,"r") as gtf,open(outdir+"/genes_loc.xls","w")as outer:
        outer.write("ChrNum\tDirector\tStarNum\tEndNum\tGeneId\tGeneFlag\tFlagValue\n")
        #loc_dict = defaultdict(set)
        for line in gtf:
            if not line.startswith("#"):
                lines = line.strip("\n").split("\t")
                chrID = lines[0]
                #sample_name = lines[0]
                flag = lines[2]
                direc = lines[6]
                star = lines[3]
                endr = lines[4]
                gene_id = re.findall(r'\sgene_id "(.*?)";',line)[0]
                try:
                    if adic and gene_id in adic.keys():
                        gene_id = adic[gene_id]
                except:
                    pass
                try:
                    motif_id = re.findall(r'\sgene_name "(.*?)";',line)[0]
                except:
                    motif_id = gene_id
                writelist = "\t".join([lines[0],lines[6],lines[3],lines[4],gene_id,lines[2],motif_id])
                outer.write(writelist+"\n")
                


def does_data_gene(data_gene):
    if len(data_gene) > 1:
        print "there is problem for the info.xls"
        exit()
    if data_gene["Director"].all()=="+":
        star = data_gene.at[0,"StarNum"]
        endr = data_gene.at[0,"EndNum"]
    else:
        endr = data_gene.at[0,"StarNum"]
        star = data_gene.at[0,"EndNum"]
    FlagValue = data_gene["FlagValue"].all()
    ChrTbl.write("\t".join([str(star),str(endr),"gene"])+"\n")
    ChrTbl.write("\t\t\tgene\t"+FlagValue+"\n")
    print "does_data_gene is over!"

def does_data_mRNA(data_mRNA):
    data_mRNA = data_mRNA.sort_values(by='StarNum')
    lenData = len(data_mRNA)
    data_mRNA.index = range(lenData)
    for n in range(lenData):
        lines = data_mRNA.loc[n]
        if lines["Director"]=="+":
            star = lines["StarNum"]
            endr = lines["EndNum"]
        else:
            endr = lines["StarNum"]
            star = lines["EndNum"]
        if n == 0:
            ChrTbl.write("\t".join([str(star),str(endr),"mRNA"])+"\n")
            FlagValue = lines["FlagValue"]
        else:
            ChrTbl.write("\t".join([str(star),str(endr)])+"\n")
    ChrTbl.write("\t\t\tproduct\t"+FlagValue+"\n")
    print "does_data_mRNA is over!"

def does_data_CDS(data_CDS):
    data_CDS = data_CDS.sort_values(by='StarNum')
    lenData = len(data_CDS)
    data_CDS.index = range(lenData)
    for n in range(lenData):
        lines = data_CDS.loc[n]
        if lines["Director"]=="+":
            star = lines["StarNum"]
            endr = lines["EndNum"]
        else:
            endr = lines["StarNum"]
            star = lines["EndNum"]
        if n == 0:
            ChrTbl.write("\t".join([str(star),str(endr),"CDS"])+"\n")
            FlagValue = lines["FlagValue"]
        else:
            ChrTbl.write("\t".join([str(star),str(endr)])+"\n")
    ChrTbl.write("\t\t\tproduct\t"+FlagValue+"\n")
    print "does_data_CDS is over!"

def does_data_exon(data_exon):
    data_exon = data_exon.sort_values(by='StarNum')
    lenData = len(data_exon)
    data_exon.index = range(lenData)
    for n in range(lenData):
        lines = data_exon.loc[n]
        if lines["Director"]=="+":
            star = lines["StarNum"]
            endr = lines["EndNum"]
        else:
            endr = lines["StarNum"]
            star = lines["EndNum"]
        ChrTbl.write("\t".join([str(star),str(endr),"exon"])+"\n")
        ChrTbl.write("\t\t\tnumber "+str(n+1)+"\n")
    print "does_data_exon is over!"

def does_data_UTR(data_UTR):
    data_UTR = data_UTR.sort_values(by='StarNum')
    lenData = len(data_UTR)
    data_UTR.index = range(lenData)
    for n in range(lenData):
        lines = data_UTR.loc[n]
        if lines["Director"]=="+":
            star = lines["StarNum"]
            endr = lines["EndNum"]
        else:
            endr = lines["StarNum"]
            star = lines["EndNum"]
        if n == 0:
            flag = lines["GeneFlag"]
            ChrTbl.write("\t".join([str(star),str(endr),flag])+"\n")
            FlagValue = lines["FlagValue"]
        else:
            ChrTbl.write("\t".join([str(star),str(endr)])+"\n")
    ChrTbl.write("\t\t\tnote\t"+FlagValue+"\n")
    print "does_data_utr is over!"

def does_data_tRNA(data_tRNA):
    data_tRNA = data_tRNA.sort_values(by='StarNum')
    lenData = len(data_tRNA)
    data_tRNA.index = range(lenData)
    for n in range(lenData):
        lines = data_tRNA.loc[n]
        if lines["Director"]=="+":
            star = lines["StarNum"]
            endr = lines["EndNum"]
        else:
            endr = lines["StarNum"]
            star = lines["EndNum"]
        flag = lines["GeneFlag"]
        if n == 0:
            ChrTbl.write("\t".join([str(star),str(endr),"tRNA"])+"\n")
            FlagValue = lines["FlagValue"]
        else:
            ChrTbl.write("\t".join([str(star),str(endr)])+"\n")
    ChrTbl.write("\t\t\tproduct\t"+FlagValue+"\n")
    print "does_data_tRNA is over!"

def does_data_ncRNA(data_ncRNA):
    data_ncRNA = data_ncRNA.sort_values(by='StarNum')
    lenData = len(data_ncRNA)
    data_ncRNA.index = range(lenData)
    for n in range(lenData):
        lines = data_ncRNA.loc[n]
        if lines["Director"]=="+":
            star = lines["StarNum"]
            endr = lines["EndNum"]
        else:
            endr = lines["StarNum"]
            star = lines["EndNum"]
        flag = lines["GeneFlag"]
        if n == 0:
            ChrTbl.write("\t".join([str(star),str(endr),flag])+"\n")
            FlagValue = lines["FlagValue"]
        else:
            ChrTbl.write("\t".join([str(star),str(endr)])+"\n")
    ChrTbl.write("\t\t\tncRNA_class\t"+FlagValue+"\n")
    ChrTbl.write("\t\t\tnote\t"+FlagValue+"\n")
    print "does_data_ncRNA is over!"

def does_data_rRNA(data_rRNA):
    azList = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
    data_rRNA = data_rRNA.sort_values(by='StarNum')
    lenData = len(data_rRNA)
    data_rRNA.index = range(lenData)
    for n in range(lenData):
        lines = data_rRNA.loc[n]
        if lines["Director"]=="+":
            star = lines["StarNum"]
            endr = lines["EndNum"]
        else:
            endr = lines["StarNum"]
            star = lines["EndNum"]
        flag = lines["GeneFlag"]
        ChrTbl.write("\t".join([str(star),str(endr),"rRNA"])+"\n")
        FlagValue = lines["FlagValue"]
        if "16S" in  FlagValue:
            ChrTbl.write("\t\t\tproduct\t16S ribosomal RNA\n")
            ChrTbl.write("\t\t\tgene\trrs{az}\n".format(az=azlist[n]))
        elif "23S" in  FlagValue:
            ChrTbl.write("\t\t\tproduct\t23S ribosomal RNA\n")
            ChrTbl.write("\t\t\tgene\trrl{az}\n".format(az=azList[n]))
        elif "5S" in  FlagValue:
            ChrTbl.write("\t\t\tproduct\t5S ribosomal RNA\n")
            ChrTbl.write("\t\t\tgene\trrf{az}\n".format(az=azList[n]))
        else:
            ChrTbl.write("\t\t\tproduct\t"+FlagValue+"\n")
            ChrTbl.write("\t\t\tgene\t"+FlagValue+"_{az}\n".format(az=azList[n]))
    print "does_data_rRNA is over!"


def does_p_data_gene(data_gene):
    if len(data_gene) > 1:
        print "there is problem for the info.xls"
        exit()
    if data_gene["Director"].all()=="+":
        star = data_gene.at[0,"StarNum"]
        endr = data_gene.at[0,"EndNum"]
    else:
        endr = data_gene.at[0,"StarNum"]
        star = data_gene.at[0,"EndNum"]
    FlagValue = data_gene["FlagValue"].all()
    ChrTbl.write("\t".join([str(star),str(endr),"gene"])+"\n")
    ChrTbl.write("\t\t\tgene\t"+FlagValue+"\n")
    ChrTbl.write("\t\t\tlocus_tag\t"+FlagValue+"\n")
    print "does_data_gene is over!"


def does_p_data_CDS(data_CDS):
    data_CDS = data_CDS.sort_values(by='StarNum')
    lenData = len(data_CDS)
    data_CDS.index = range(lenData)
    for n in range(lenData):
        lines = data_CDS.loc[n]
        if lines["Director"]=="+":
            star = lines["StarNum"]
            endr = lines["EndNum"]
        else:
            endr = lines["StarNum"]
            star = lines["EndNum"]
        if n == 0:
            ChrTbl.write("\t".join([str(star),str(endr),"CDS"])+"\n")
            FlagValue = lines["FlagValue"]
        else:
            ChrTbl.write("\t".join([str(star),str(endr)])+"\n")
    if args.lab:
        ChrTbl.write("\t\t\tproduct\tgnl|"+args.lab+"|"+FlagValue+"\n")
    else:
        pass
    ChrTbl.write("\t\t\tproduct_id\t"+FlagValue+"\n")
    print "does_data_CDS is over!"


def does_p_data_tRNA(data_tRNA):
    data_tRNA = data_tRNA.sort_values(by='StarNum')
    lenData = len(data_tRNA)
    data_tRNA.index = range(lenData)
    for n in range(lenData):
        lines = data_tRNA.loc[n]
        if lines["Director"]=="+":
            star = lines["StarNum"]
            endr = lines["EndNum"]
        else:
            endr = lines["StarNum"]
            star = lines["EndNum"]
        flag = lines["GeneFlag"]
        if n == 0:
            ChrTbl.write("\t".join([str(star),str(endr),"tRNA"])+"\n")
            FlagValue = lines["FlagValue"]
            GeneId = lines["GeneId"]
        else:
            ChrTbl.write("\t".join([str(star),str(endr)])+"\n")
    ChrTbl.write("\t\t\tproduct\t"+FlagValue+"\n")
    ChrTbl.write("\t\t\tgene\t"+GeneId+"\n")
    print "does_data_tRNA is over!"

def does_p_data_ncRNA(data_ncRNA):
    data_ncRNA = data_ncRNA.sort_values(by='StarNum')
    lenData = len(data_ncRNA)
    data_ncRNA.index = range(lenData)
    for n in range(lenData):
        lines = data_ncRNA.loc[n]
        if lines["Director"]=="+":
            star = lines["StarNum"]
            endr = lines["EndNum"]
        else:
            endr = lines["StarNum"]
            star = lines["EndNum"]
        flag = lines["GeneFlag"]
        if n == 0:
            ChrTbl.write("\t".join([str(star),str(endr),flag])+"\n")
            FlagValue = lines["FlagValue"]
            GeneId = lines["GeneId"]
        else:
            ChrTbl.write("\t".join([str(star),str(endr)])+"\n")
    ChrTbl.write("\t\t\tncRNA_class\t"+FlagValue+"\n")
    ChrTbl.write("\t\t\tnote\t"+FlagValue+"\n")
    ChrTbl.write("\t\t\tgene\t"+GeneId+"\n")
    print "does_data_ncRNA is over!"

def geneanno(cfile):
    #print cfile
    fr=open(cfile)
    adic={}
    for i in fr:
	i=i.strip()
	spi=i.split('\t')
	if spi[0] not in adic.keys():
	    adic[spi[0]]=spi[1]
    return adic

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="huoqu can shu")
    parser.add_argument('-o', '--out',required=False, default = './',help='ouput dir ')
    parser.add_argument('-s', '--tem',required=True, help='the file of template.sbt')
    parser.add_argument('-f', '--fasta',required=True, help='the file of fasta.fsa')
    parser.add_argument('-t', '--tbl',required=False, help='the tbl file ')
    parser.add_argument('-i', '--info',required=False, help=' file info')
    parser.add_argument('-c', '--clas',required=False,default='Eukaryotes', help=' class for your specise:Eukaryotes or Prokaryotes .default == Eukaryotes')
    parser.add_argument('-l', '--lab',required=False, help=' the name of your laboratory')
    parser.add_argument('-g', '--gtf',required=False, help=' file with gtf')
    parser.add_argument('-a', '--anno',required=False, help='the annotation file for moti_name')
    parser.add_argument('--tbl2asn',required=True, help='tbl2asn path')
    args = parser.parse_args()
    tbl2asn = args.tbl2asn
    temsbt = os.path.abspath(args.tem)
    fasta = os.path.abspath(args.fasta)
    outdir = os.path.abspath(args.out)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    if not args.tbl:
        if not args.info:
            if not args.gtf:
                print "lack of file with xxx.tbl need"
                jobs = "{tbl2asn} -i {fasta} -o {outdir}/test.sqn  -a s".format(tbl2asn=tbl2asn,fasta=fasta,outdir=outdir)
                #subprocess.call(jobs,shell=True)
                os.system(jobs)
            else:
                if args.anno:
                    adic=geneanno(args.anno)
                gtfile = os.path.abspath(args.gtf)
                arrange_gtf(gtfile,outdir)
                args.info = outdir+"/genes_loc.xls"
        #else:
        dataAll = pd.read_csv(args.info,header=0,sep="\t")
        chrList = list(set(dataAll["ChrNum"]))
        fristChr = chrList[0]
        for eachChr in chrList:
            #print eachChr
            ChrTbl = open("{outdir}/{chr}.tbl".format(outdir=outdir,chr=eachChr),"w")
            ChrTbl.write(">Feature {chr}\n".format(chr=eachChr))
            dataChr = dataAll.loc[np.where(dataAll["ChrNum"]==eachChr)]
            dataChr.index = range(len(dataChr))
            #print dataChr.head()
            geneList = list(set(dataChr["GeneId"]))
            if args.clas == "Eukaryotes":
                for eachGene in geneList:
                    dataGene = dataChr.loc[np.where(dataChr["GeneId"]==eachGene)]
                    dataGene.index = range(len(dataGene))
                    try:
                        data_gene = dataGene.loc[np.where(dataGene["GeneFlag"]=="gene")]
                        does_data_gene(data_gene)
                    except:
                        print "There is a problem for the {info}".format(info=args.info)
                    try:
                        #print "mRNA"
                        data_mRNA = dataGene.loc[np.where(dataGene["GeneFlag"]=="mRNA")]
                        does_data_mRNA(data_mRNA)
                    except:
                        #print "wrong"
                        pass
                    try:
                        #print "CDS"
                        data_CDS = dataGene.loc[np.where(dataGene["GeneFlag"]=="CDS")]
                        does_data_CDS(data_CDS)
                    except:
                        #print "wrong"
                        pass
                    try:
                        #print "exon"
                        data_exon = dataGene.loc[np.where(dataGene["GeneFlag"]=="exon")]
                        does_data_exon(data_exon)
                    except:
                        #print "wrong"
                        pass
                    try:
                        #print "UTR"
                        data_UTR = dataGene.loc[np.where((dataGene["GeneFlag"]=="UTR")|(dataGene["GeneFlag"]=="5'UTR")|(dataGene["GeneFlag"]=="3'UTR"))]
                        does_data_UTR(data_UTR)
                    except:
                        #print "wrong"
                        pass
                    try:
                        #print "tRNA"
                        data_tRNA = dataGene.loc[np.where(dataGene["GeneFlag"]=="tRNA")]
                        does_data_tRNA(data_tRNA)
                    except:
                        #print "wrong"
                        pass
                    try:
                        #print "ncRNA"
                        data_ncRNA = dataGene.loc[np.where(dataGene["GeneFlag"]=="ncRNA")]
                        does_data_ncRNA(data_ncRNA)
                    except:
                        #print "wrong"
                        pass
                    try:
                        #print "rRNA"
                        data_rRNA = dataGene.loc[np.where(dataGene["GeneFlag"]=="rRNA")]
                        does_data_rRNA(data_rRNA)
                    except:
                        #print "wrong"
                        pass
            elif args.clas == "Prokaryotes":
                for eachGene in geneList:
                    dataGene = dataChr.loc[np.where(dataChr["GeneId"]==eachGene)]
                    dataGene.index = range(len(dataGene))
                    try:
                        data_gene = dataGene.loc[np.where(dataGene["GeneFlag"]=="gene")]
                        does_p_data_gene(data_gene)
                    except:
                        print "There is a problem for the {info}".format(info=args.info)
                    try:
                        #print "CDS"
                        data_CDS = dataGene.loc[np.where(dataGene["GeneFlag"]=="CDS")]
                        if len(data_CDS) == 0:
                            print "Prokaryotes could not lack CDS"
                            exit()
                        does_p_data_CDS(data_CDS)
                    except:
                        #print "wrong"
                        pass
                    try:
                        #print "UTR"
                        data_UTR = dataGene.loc[np.where((dataGene["GeneFlag"]=="UTR")|(dataGene["GeneFlag"]=="5'UTR")|(dataGene["GeneFlag"]=="3'UTR"))]
                        does_data_UTR(data_UTR)
                    except:
                        #print "wrong"
                        pass
                    try:
                        #print "tRNA"
                        data_tRNA = dataGene.loc[np.where(dataGene["GeneFlag"]=="tRNA")]
                        does_p_data_tRNA(data_tRNA)
                    except:
                        #print "wrong"
                        pass
                    try:
                        #print "ncRNA"
                        data_ncRNA = dataGene.loc[np.where((dataGene["GeneFlag"]=="ncRNA")|(dataGene["GeneFlag"]=="misc_RNA"))]
                        if len(data_ncRNA) == 0:
                            print "Prokaryotes could not lack ncRNA"
                            exit()
                        does_p_data_ncRNA(data_ncRNA)
                    except:
                        #print "wrong"
                        pass
                    try:
                        #print "rRNA"
                        data_rRNA = dataGene.loc[np.where(dataGene["GeneFlag"]=="rRNA")]
                        if len(data_rRNA) == 0:
                            print "Prokaryotes could not lack rRNA"
                            exit()
                        does_p_data_rRNA(data_rRNA)
                    except:
                        #print "wrong"
                        pass
            ChrTbl.close()
    else:
        tblFile = os.path.abspath(args.tbl)
        fristChr  = tblFile.split("/")[-1]
        print "#####1#####"
        #subprocess.call("cp {tblFile} {outdir}/{fristChr}.tbl".format(tblFile=tblFile,outdir=outdir,fristChr=fristChr),shell=True)
        os.system("cp {tblFile} {outdir}/{fristChr}.tbl".format(tblFile=tblFile,outdir=outdir,fristChr=fristChr))
    print "#####2#####"
    fasjob = "cp {fasta} {outdir}/{fristChr}.fsa".format(fasta=fasta,fristChr=fristChr,outdir=outdir)
    #subprocess.call(fasjob,shell=True)
    os.system(fasjob)
    tbljobs = "{tbl2asn}  -t {template}  -p {outdir} -a s -V vb".format(tbl2asn=tbl2asn,template=temsbt,outdir=outdir)
    print "#####3#####"
    #subprocess.call(tbljobs,shell=True)
    os.system(tbljobs)
    print "#####4#####"
