# -*- coding: utf-8 -*-
from collections import defaultdict
from collections import OrderedDict
import os
import os
import re
import logging
import numpy as np
import pandas as pd
import argparse
import subprocess


class GffFormat(object):
    def __init__(self):
        self.genelist=[]
        self.translist=[]
        self.gene2trans=defaultdict(list)
        # self.gene2name=defaultdict(str)
        # self.trans2gene=defaultdict(str)
        # self.trans2name = defaultdict(str)
        self.gene2name = dict()
        self.trans2gene = dict()
        self.trans2name = dict()
        self.id2pos = defaultdict(list)
        self.name2pos = defaultdict(list)

    def get_info(self,gff=None,gtf=None):
        with open(gff,"r") as g,open(gtf,"w") as t:
            for line in g.readlines():
              if not line.startswith("#"):
                  seqid, source, type, start, end, score, strand, phase, order_detail_attributes =self.detail_analysis(line)
                  t.write("\t".join([seqid, source, type, start, end, score, strand, phase])+"\t")
                  for key in order_detail_attributes:
                      t.write(" "+key+" \""+order_detail_attributes[key]+"\";")
                  t.write("\n")

    def detail_analysis(self,line):
             term_infos = line.strip().split("\t")
             seqid, source, type, start, end, score, strand, phase, attributes = \
                 term_infos[0], term_infos[1], term_infos[2], term_infos[3], term_infos[4], term_infos[5], term_infos[6], term_infos[7], term_infos[8]
             if "gene" in type:
                 g=re.search(r'gene_id=(.+?);', attributes)
                 if g:
                     gene_id = g.group(1)
                 elif re.search(r'ID=gene:(.+?);', attributes):
                     gene_id=re.search(r'ID=gene:(.+?);', attributes).group(1)
                 else:
                     transcript_id=re.search(r'ID=transcript:(.+?);', attributes).group(1)
                     gene_id = re.search(r'Parent=gene:(.+?);', attributes).group(1)
                 if  not 'transcript_id' in locals().keys():
                    self.genelist.append(gene_id)
                    gene_name = re.search(r'Name=(.+?);', attributes).group(1)
                    self.gene2name[gene_id] = gene_name
                 else:
                     transcript_name=re.search(r'Name=(.+?);', attributes).group(1)
                     try:
                         self.trans2name[transcript_id] = transcript_name
                         gene_name= self.gene2name[gene_id]
                     except:
                         pass
                     try:
                         self.trans2gene[transcript_id] = gene_id
                     except:
                         pass
                 detail_attributes = {x: y for x, y in [x.split("=") for x in attributes.split(";")]}
                 try:
                     del detail_attributes["ID"]
                 except:
                     pass
                 try:
                     del detail_attributes["Name"]
                 except:
                     pass
                 try:
                     del detail_attributes["gene_id"]
                 except:
                     pass
                 order_detail_attributes = OrderedDict()
                 try:
                    order_detail_attributes["gene_id"] = gene_id
                 except:
                     pass
                 if 'transcript_id' in locals().keys():
                     order_detail_attributes["transcript_id"] = transcript_id
                     order_detail_attributes["transcript_name"] = transcript_name
                 try:
                    order_detail_attributes["gene_name"] = gene_name
                 except:
                     pass
                 for i in detail_attributes:
                     order_detail_attributes[i] = detail_attributes[i]

             elif type == "lnc_RNA" or  type == "mRNA" or  type == "transcript" or "transcript" in type:
                 term_infos = line.strip().split("\t")
                 t=re.search(r'transcript_id=(.+?);', attributes)
                 try:
                     transcript_id = t.group(1) if t else re.search(r'ID=transcript:(.+?);', attributes).group(1)
                 except:
                     pass
                 try:
                    gene_id=re.search(r'Parent=gene:(.+?);', attributes).group(1)
                    self.trans2gene[transcript_id] = gene_id
                 except:
                     pass
                 try:
                    gene_name=self.gene2name[gene_id]
                 except:
                     logging.info("{}没有gene_name".format(gene_id))
                 transcript_name=re.search(r'Name=(.+?);', attributes).group(1)
                 try:
                    self.trans2name[transcript_id] = transcript_name
                 except:
                     pass
                 detail_attributes = {x: y for x, y in [x.split("=") for x in attributes.split(";")]}
                 try:
                     del detail_attributes["ID"]
                 except:
                     pass
                 try:
                     del detail_attributes["Parent"]
                 except:
                     pass
                 try:
                     del detail_attributes["Name"]
                 except:
                     pass
                 try:
                     del detail_attributes["transcript_id"]
                 except:
                     pass
                 order_detail_attributes = OrderedDict()
                 try:
                     order_detail_attributes["gene_id"] = gene_id
                 except:
                     pass
                 try:
                    order_detail_attributes["transcript_id"] = transcript_id
                 except:
                     pass
                 try:
                     order_detail_attributes["gene_name"] = gene_name
                 except:
                     pass
                 try:
                     order_detail_attributes["transcript_name"] = transcript_name
                 except:
                     pass
                 for i in detail_attributes:
                     order_detail_attributes[i] = detail_attributes[i]

             elif type == "exon":
                 term_infos = line.strip().split("\t")
                 e=re.search(r'exon_id=(.+?);', attributes)
                 try:
                    exo_id = e.group(1) if e else re.search(r'ID=id:(.+?);', attributes).group(1)
                 except:
                     pass
                 try:
                    transcript_id = re.search(r'Parent=transcript:(.+?);', attributes).group(1)
                 except:
                     pass
                 try:
                    gene_id = self.trans2gene[transcript_id]
                 except:
                     pass
                 try:
                    gene_name = self.gene2name[gene_id]
                 except:
                     pass
                 try:
                    transcript_name=self.trans2name[transcript_id]
                 except:
                     pass
                 detail_attributes = {x: y for x, y in [x.split("=") for x in attributes.split(";")]}
                 try:
                     del detail_attributes["Parent"]
                 except:
                     pass
                 try:
                     del detail_attributes["ID"]
                 except:
                     pass
                 try:
                     del detail_attributes["Name"]
                 except:
                     pass
                 try:
                     del detail_attributes["exon_id"]
                 except:
                     pass
                 order_detail_attributes = OrderedDict()
                 try:
                    order_detail_attributes["gene_id"] = gene_id
                 except:
                     pass
                 try:
                    order_detail_attributes["transcript_id"] = transcript_id
                 except:
                     pass
                 try:
                    order_detail_attributes["exo_id"] = exo_id
                 except:
                        pass
                 try:
                     order_detail_attributes["gene_name"] = gene_name
                 except:
                        pass
                 try:
                     order_detail_attributes["transcript_name"] = transcript_name
                 except:
                        pass
                 for i in detail_attributes:
                     order_detail_attributes[i] = detail_attributes[i]


             elif type == "CDS":
                 term_infos = line.strip().split("\t")
                 e = re.search(r'protein_id=(.+?);', attributes)
                 try:
                    protein_id = e.group(1) if e else re.search(r'ID=CDS:(.+?);', attributes).group(1)
                 except:
                        pass
                 try:
                    transcript_id = re.search(r'Parent=transcript:(.+?);', attributes).group(1)
                 except:
                        pass
                 try:
                    gene_id = self.trans2gene[transcript_id]
                 except:
                        pass
                 try:
                    gene_name = self.gene2name[gene_id]
                 except:
                        pass
                 try:
                    transcript_name = self.trans2name[transcript_id]
                 except:
                        pass
                 detail_attributes = {x: y for x, y in [x.split("=") for x in attributes.split(";")]}
                 try:
                     del detail_attributes["Parent"]
                 except:
                     pass
                 try:
                     del detail_attributes["ID"]
                 except:
                     pass
                 order_detail_attributes = OrderedDict()
                 try:
                    order_detail_attributes["gene_id"] = gene_id
                 except:
                        pass
                 try:
                    order_detail_attributes["transcript_id"] = transcript_id
                 except:
                        pass
                 try:
                     order_detail_attributes["gene_name"] = gene_name
                 except:
                        pass
                 try:
                     order_detail_attributes["transcript_name"] = transcript_name
                 except:
                        pass
                 try:
                     order_detail_attributes["protein_id"] = protein_id
                 except:
                        pass
                 for i in detail_attributes:
                     order_detail_attributes[i] = detail_attributes[i]

             elif type == "rRNA":
                 term_infos = line.strip().split("\t")
                 try:
                    t = re.search(r'transcript_id=(.+?);', attributes).group(1)
                 except:
                        pass
                 try:
                    transcript_id = t if t else re.search(r'ID=transcript:(.+?);', attributes).group(1)
                 except:
                        pass
                 try:
                    gene_id = re.search(r'Parent=gene:(.+?);', attributes).group(1)
                 except:
                        pass
                 try:
                    self.trans2gene[transcript_id] = gene_id
                 except:
                        pass
                 try:
                    gene_name = self.gene2name[gene_id]
                 except:
                        pass
                 try:
                    transcript_name =re.search(r'Name=(.+?);', attributes).group(1) if re.search(r'Name=(.+?);', attributes).group(1)\
                        else self.trans2name[transcript_id]
                 except:
                        pass
                 detail_attributes = {x: y for x, y in [x.split("=") for x in attributes.split(";")]}
                 try:
                     del detail_attributes["Parent"]
                 except:
                     pass
                 try:
                     del detail_attributes["ID"]
                 except:
                     pass
                 try:
                     del detail_attributes["Name"]
                 except:
                     pass
                 try:
                     del detail_attributes["transcript_id"]
                 except:
                     pass
                 order_detail_attributes = OrderedDict()
                 try:
                    order_detail_attributes["gene_id"] = gene_id
                 except:
                        pass
                 try:
                    order_detail_attributes["transcript_id"] = transcript_id
                 except:
                        pass
                 try:
                     order_detail_attributes["gene_name"] = gene_name
                 except:
                        pass
                 try:
                     order_detail_attributes["transcript_name"] = transcript_name
                 except:
                        pass
                 for i in detail_attributes:
                     order_detail_attributes[i] = detail_attributes[i]


             elif type == "tRNA":
                 term_infos = line.strip().split("\t")
                 t = re.search(r'transcript_id=(.+?);', attributes)
                 try:
                     transcript_id = t.group(1) if t else re.search(r'ID=transcript:(.+?);', attributes).group(1)
                 except:
                        pass
                 try:
                     gene_id = re.search(r'Parent=gene:(.+?);', attributes).group(1)
                 except:
                        pass
                 try:
                    self.trans2gene[transcript_id] = gene_id
                 except:
                        pass
                 try:
                    gene_name = self.gene2name[gene_id]
                 except:
                        pass
                 try:
                     transcript_name =re.search(r'Name=(.+?);', attributes).group(1) if re.search(r'Name=(.+?);', attributes).group(1)\
                         else self.trans2name[transcript_id]
                 except:
                        pass
                 detail_attributes = {x: y for x, y in [x.split("=") for x in attributes.split(";")]}
                 try:
                     del detail_attributes["Parent"]
                 except:
                     pass
                 try:
                     del detail_attributes["ID"]
                 except:
                     pass
                 try:
                     del detail_attributes["Name"]
                 except:
                     pass
                 try:
                     del detail_attributes["transcript_id"]
                 except:
                     pass
                 order_detail_attributes = OrderedDict()
                 try:
                    order_detail_attributes["gene_id"] = gene_id
                 except:
                        pass
                 try:
                    order_detail_attributes["transcript_id"] = transcript_id
                 except:
                        pass
                 try:
                     order_detail_attributes["gene_name"] = gene_name
                 except:
                        pass
                 try:
                     order_detail_attributes["transcript_name"] = transcript_name
                 except:
                        pass
                 for i in detail_attributes:
                     order_detail_attributes[i] = detail_attributes[i]

             elif "RNA" in type:
                 term_infos = line.strip().split("\t")
                 try:
                    t = re.search(r'transcript_id=(.+?);', attributes).group(1)
                 except:
                        pass
                 try:
                    transcript_id = t if t else re.search(r'ID=transcript:(.+?);', attributes).group(1)
                 except:
                        pass
                 try:
                    gene_id = re.search(r'Parent=gene:(.+?);', attributes).group(1)
                 except:
                        pass
                 try:
                    self.trans2gene[transcript_id] = gene_id
                 except:
                        pass
                 try:
                    gene_name = self.gene2name[gene_id]
                 except:
                        pass
                 try:
                     transcript_name =re.search(r'Name=(.+?);', attributes).group(1) if re.search(r'Name=(.+?);', attributes).group(1)\
                         else self.trans2name[transcript_id]
                 except:
                        pass
                 detail_attributes = {x: y for x, y in [x.split("=") for x in attributes.split(";")]}
                 try:
                     del detail_attributes["Parent"]
                 except:
                     pass
                 try:
                     del detail_attributes["ID"]
                 except:
                     pass
                 try:
                     del detail_attributes["Name"]
                 except:
                     pass
                 try:
                     del detail_attributes["transcript_id"]
                 except:
                     pass
                 order_detail_attributes = OrderedDict()
                 try:
                    order_detail_attributes["gene_id"] = gene_id
                 except:
                        pass
                 try:
                    order_detail_attributes["transcript_id"] = transcript_id
                 except:
                        pass
                 try:
                     order_detail_attributes["gene_name"] = gene_name
                 except:
                        pass
                 try:
                     order_detail_attributes["transcript_name"] = transcript_name
                 except:
                        pass
                 for i in detail_attributes:
                     order_detail_attributes[i] = detail_attributes[i]


             else:
                 term_infos = line.strip().split("\t")
                 t = re.search(r'Parent=transcript:(.+?);', attributes)
                 if t:
                     transcript_id = t.group(1)
                     try:
                         transcript_name =  self.trans2name[transcript_id]
                     except:
                         pass
                     try:
                        gene_id = re.search(r'Parent=gene:(.+?);', attributes).group(1)
                     except:
                         pass
                     try:
                        self.trans2gene[transcript_id] = gene_id
                     except:
                         pass
                     try:
                        gene_name = self.gene2name[gene_id]
                     except:
                         pass
                     order_detail_attributes = OrderedDict()
                     try:
                        order_detail_attributes["gene_id"] = gene_id
                     except:
                         pass
                     try:
                        order_detail_attributes["transcript_id"] = transcript_id
                     except:
                         pass
                     try:
                         order_detail_attributes["gene_name"] = gene_name
                     except:
                         pass
                     try:
                         order_detail_attributes["transcript_name"] = transcript_name
                     except:
                         pass
                 else:
                     order_detail_attributes=OrderedDict()


             return seqid, source, type, start, end, score, strand, phase, order_detail_attributes

    def bed_by_id_list(self, gff,list,bed):

        with open(gff,"r") as g:
            for line in g.readlines():
                if not line.startswith("#"):
                    line=line.strip().split("\t")
                    chorm,type, start, end, score, strand, phase, attributes = line[0],line[2],line[3],line[4],line[5],line[6],line[7],line[8]
                    if type == "gene":
                        gene_id = re.search(r'ID=gene:(.+?);', attributes).group(1)
                        gene_name = re.search(r'Name=(.+?);', attributes).group(1)
                        self.id2pos[gene_id]=[chorm,start,end,score,strand,phase]
                        self.name2pos[gene_name] = [chorm,start, end,score, strand, phase]
        with open(list,"r") as l,open(bed ,"w") as b,open("no_exist_list","w") as n:
            for line in l.readlines():
                name = line.strip()
                if name in self.id2pos:
                    detail=[name]
                    detail.extend(self.id2pos[name])
                    b.write("\t".join(detail)+"\n")
                elif name in self.name2pos:
                    detail = [name]
                    detail.extend(self.name2pos[name])
                    b.write("\t".join(detail) + "\n")
                else:
                    n.write(name+"\n")

    def gfftobed(self,gff,bed):

        jobs = "convert2bed --input=gff --output=bed < {} > {}".format(gff,bed)
        subprocess.call(jobs, shell=True)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-gff', type=str, required=True,
                        help='gff file for analysis.')
    parser.add_argument('-genelist', type=str, default=None,
                        help="genelist to extract positions from fasta  file"
                             )
    parser.add_argument('-out', type=str, required=True,help="output_dir")
    parser.add_argument('-type', type=str, required=True, help="processing algorithm")
    args = parser.parse_args()
    gff_file = args.gff
    gff_name=os.path.splitext(os.path.basename(gff_file))[0]
    type = args.type
    output_dir = args.out
    gff = GffFormat()
    if type == "gfftogtf":
        gtf_file=os.path.join(output_dir,gff_name+".gtf")
        gff.get_info(gff=gff_file,gtf= gtf_file)
    elif type == "extract_pos":
        gene_list = args.genelist
        bed_path=os.path.join(output_dir,"extrac.bed")
        gff.bed_by_id_list(gff_file,gene_list,bed_path)
    elif type == "gfftobed":
        outbed=os.path.join(output_dir,"extrac.bed")
        gff.gfftobed(gff_file,outbed)

    # gff.get_info(gff="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab_batch2/gtf_bed_gff/data/gff/Homo_sapiens.GRCh38.99.gff3",
    #              gtf="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab_batch2/gtf_bed_gff/data/gff/testself_new.gtf")


    # gff.bed_by_id_list("/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab_batch2/gtf_bed_gff/data/gff/Homo_sapiens.GRCh38.99.gff3",
    #                    "/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab_batch2/gtf_bed_gff/data/gff/testlist",
    #                    "/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab_batch2/gtf_bed_gff/data/gff/extrac.bed")











