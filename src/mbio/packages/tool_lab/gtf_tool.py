# -*- coding: utf-8 -*-
from collections import defaultdict
from collections import OrderedDict
import os
import os
import re
import logging
import argparse
import subprocess

class GtfFormat(object):
    def __init__(self):
        self.genelist=[]
        self.translist=[]
        self.gene2trans=defaultdict(list)
        self.gene2name=defaultdict(str)
        self.trans2gene=defaultdict(str)
        self.trans2name = defaultdict(str)

    def gtf_to_gff(self,gtf=None,gff=None):
        with open(gtf,"r") as g,open(gff,"w") as t:
            for line in g.readlines():
              if not line.startswith("#"):
                  seqid, source, type, start, end, score, strand, phase, order_detail_attributes =self.detail_analysis(line)
                  t.write("\t".join([seqid, source, type, start, end, score, strand, phase])+"\t")
                  for key in order_detail_attributes:
                      t.write(key+"="+order_detail_attributes[key]+";")
                  t.write("\n")

    def detail_analysis(self,line):
             term_infos = line.strip().split("\t")
             seqid, source, type, start, end, score, strand, phase, attributes = \
                 term_infos[0], term_infos[1], term_infos[2], term_infos[3], term_infos[4], term_infos[5], term_infos[6], term_infos[7], term_infos[8]
             if "gene" in type:
                 g = re.search(r'gene_id \"(.+?)\";', attributes)
                 if g:
                     gene_id = g.group(1)
                     self.genelist.append(gene_id)
                     # detail_attributes = {x.strip() : y.strip("\"") for x, y in [x.strip().split(" ") for x in attributes.strip(";").split(";")]}
                     detail_attributes = {x.strip(): y.strip("\"") for x, y in
                                          [x.strip().split(" \"") for x in attributes.strip(";").split(";")]}
                     try:
                         del detail_attributes["gene_id"]
                     except:
                         pass
                     order_detail_attributes = OrderedDict()
                     order_detail_attributes["ID"] = gene_id
                     for i in detail_attributes:
                         order_detail_attributes[i] = detail_attributes[i]
                 else:
                     order_detail_attributes = OrderedDict()


             elif type == "transcript":
                 t=re.search(r'transcript_id \"(.+?)\";', attributes)
                 if t:
                     transcript_id = t.group(1)
                     g = re.search(r'gene_id \"(.+?)\";', attributes)
                     if g:
                         gene_id = g.group(1)
                     try:
                        self.trans2gene[transcript_id] = gene_id
                     except:
                         pass
                     detail_attributes = {x.strip(): y.strip("\"") for x, y in
                                          [x.strip().split(" \"") for x in attributes.strip(";").split(";")]}
                     # detail_attributes = {x.strip(): y.strip("\"") for x,y in
                     #                      [x.strip().split(" \"") for x in attributes.strip(";").split(";")]}
                     order_detail_attributes = OrderedDict()
                     try:
                        order_detail_attributes["ID"]=transcript_id
                     except:
                         pass
                     try:
                        order_detail_attributes["Parent"] = gene_id
                     except:
                         pass
                     for i in detail_attributes:
                         order_detail_attributes[i] = detail_attributes[i]

                 else:
                     order_detail_attributes = OrderedDict()

             elif type == "CDS" or type == "exon":
                 t = re.search(r'transcript_id \"(.+?)\";', attributes)
                 if t:
                     transcript_id = t.group(1)
                     # detail_attributes = {x.strip(): y.strip("\"") for x, y in
                     #                      [x.strip().split(" ") for x in attributes.strip(";").split(";")]}
                     detail_attributes = {x.strip(): y.strip("\"") for x, y in
                                          [x.strip().split(" \"") for x in attributes.strip(";").split(";")]}
                     order_detail_attributes = OrderedDict()
                     order_detail_attributes["Parent"] = transcript_id
                     for i in detail_attributes:
                         order_detail_attributes[i] = detail_attributes[i]
                 else:
                     order_detail_attributes = OrderedDict()
                     detail_attributes = {x.strip(): y.strip("\"") for x, y in
                                          [x.strip().split(" \"") for x in attributes.strip(";").split(";")]}
                     for i in detail_attributes:
                         order_detail_attributes[i] = detail_attributes[i]


             else:
                 t = re.search(r'transcript_id \"(.+?)\";', attributes)
                 if t:
                     transcript_id = t.group(1)
                     # detail_attributes = {x.strip(): y.strip(" \"") for x, y in
                     #                      [x.strip().split(" ") for x in attributes.strip(";").split(";")]}
                     detail_attributes = {x.strip(): y.strip("\"") for x, y in
                                          [x.strip().split(" \"") for x in attributes.strip(";").split(";")]}
                     order_detail_attributes = OrderedDict()
                     order_detail_attributes["Parent"] = transcript_id
                     for i in detail_attributes:
                         order_detail_attributes[i] = detail_attributes[i]
                 else:
                     order_detail_attributes = OrderedDict()

             return seqid, source, type, start, end, score, strand, phase, order_detail_attributes

    def gtftobed(self, gtf, bed):
        jobs = "convert2bed --input=gff --output=bed < {} > {}".format(gtf, bed)
        subprocess.call(jobs, shell=True)




if __name__ == "__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument('-gtf', type=str, required=True,
                            help='gtf file for analysis.')
        # parser.add_argument('-genelist', type=str, default=None,
        #                     help="genelist to extract positions from fasta  file"
        #                     )
        parser.add_argument('-out', type=str, required=True, help="output_dir")
        parser.add_argument('-type', type=str, required=True, help="processing algorithm")
        args = parser.parse_args()
        gtf_file = args.gtf
        gtf_name = os.path.splitext(os.path.basename(gtf_file))[0]
        type = args.type
        output_dir = args.out
        gtf = GtfFormat()
        if type == "gtftogff":
            gff_file = os.path.join(output_dir, gtf_name+".gff")
            gtf.gtf_to_gff(gtf=gtf_file, gff=gff_file)
        elif type == "gtftobed":
            outbed = os.path.join(output_dir, "extract.bed")
            gtf.gtftobed(gtf_file, outbed)
        # elif type == "extract_pos":
        #     gene_list = args.genelist
        #     bed_path = os.path.join(output_dir, "extrac.bed")
        #     gff.bed_by_id_list(gff_file, gene_list, bed_path)

    # gff = GtfFormat()
    # gff.gtf_to_gff(gtf="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab_batch2/gtf_bed_gff/data/gtf/Homo_sapiens.GRCh38.99.gtf",
    #              gff="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab_batch2/gtf_bed_gff/data/gtf/test.gff")







