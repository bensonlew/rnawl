# -*- coding: utf-8 -*-
import os
import math
import numpy as np
import pandas as pd
import argparse
import glob
import re
import gzip
import json
import copy
import csv
import pandas as pd
from collections import defaultdict
from collections import OrderedDict
from mbio.api.database.gene_db import genome



class MergeGene(object):
    def __init__(self, genome_acc):
        # self.genome_api = genome.Genome(None)
        self.genome_acc = genome_acc
        self.gene2std = dict()
        self.std2gene = dict()
        self.use_enterz_genes = list()
        self.refgene2stdgene = dict()
        self.stdgene_dict = dict()
        # self.ref2ens_keys = OrderedDict({})

        self.ref2ens_keys = OrderedDict({
            "entrezgene_id": "GeneID"
        })
        self.ref2ens_t_keys = OrderedDict({
            "refseq_peptide": 'protein_id',
            "refseq_ncrna": 'transcript_id',
            "refseq_peptide_predicted": 'protein_id',
            "refseq_ncrna_predicted": 'transcript_id',
            "refseq_mrna": 'transcript_id',
            "refseq_dna": 'transcript_id',
            "refseq_mrna_predicted": 'transcript_id'
        })

        self.use_refseq_fields = {"entrezgene_id": "GeneID",
                                  "entrezgene_accession": "Name",
                                  "entrezgene_description": "description",
                                  "gene_biotype": "gene_biotype"}

        self.use_refseq_t_fields = {"entrezgene_id": "GeneID",
                                    "entrezgene_accession": "Name",
                                    "entrezgene_description": "description",
                                    "gene_biotype": "gene_biotype",
                                    "refseq_peptide": 'protein_id',
                                    "refseq_mrna": 'transcript_id',
                                    "refseq_ncrna": 'transcript_id'}

    def merge_pipe(self, gene_table, tran_table, ensemble_gene, ensemble_tran):
        # chr_list = self.get_genome_info(self.genome_acc)
        gene_dict, gene_id2std = self.get_gene_info(gene_table)
        # self.std2enterz = gene_id2std
        self.filter_ensemble_gene(gene_dict, gene_id2std,
                                  ensemble_gene, ensemble_gene + ".merge_refseq")

        tran_dict, tran_id2std = self.get_tran_info(tran_table)
        self.filter_ensemble_tran(tran_dict, tran_id2std,
                                   ensemble_tran, ensemble_tran + ".merge_refseq")

    def get_tran_info(self, tran_table=None):
        tran_dict = OrderedDict()
        tran_id2std = dict()
        tran_std_id = dict()

        with open(tran_table, "r") as f:
            for dic in csv.DictReader(f, delimiter='\t'):
                if dic["sg_gene_id"] in self.refgene2stdgene:
                    dic["sg_gene_id"] = self.refgene2stdgene[dic["sg_gene_id"]]
                if dic["transcript_id"] != "":
                    dic["transcript_id"] = dic["transcript_id"].split(".")[0]
                if dic["protein_id"] != "":
                    dic["protein_id"] = dic["protein_id"].split(".")[0]
                if dic["sg_gene_id"] in tran_dict:
                    tran_dict[dic["sg_gene_id"]].update({dic["sg_tran_id"]: dic})
                else:
                    tran_dict[dic["sg_gene_id"]] ={dic["sg_tran_id"]: dic}
        return tran_dict, tran_id2std

    def get_gene_info(self, gene_table=None, feature_table=None):
        gene_dict = OrderedDict()
        gene_id2std = dict()
        gene_std_id = dict()

        with open(gene_table, "r") as f:
            for dic in csv.DictReader(f, delimiter='\t'):
                gene_dict[dic["sg_gene_id"]] = dic
                gene_id2std[dic["GeneID"]] = dic["sg_gene_id"]
        return gene_dict, gene_id2std

        '''
        # 由 feature_table 获取所需信息
        if feature_table:
            with gzip.open(feature_table, 'r') as f:
                header = f.readline()
                name_old = header.lstrip("# ").split("\t")
                name_std = [re.sub(r'[^\da-zA-Z_]', "_", name).lower() for name in name_old]
                print dict(zip(name_old, name_std))

                for line in f:
                    cols = line.strip("\n").split("\t")
                    record_dict = dict(zip(name_std, cols))
                    record_dict['product_accession'] = record_dict['product_accession'].split(".")[0]
                    record_dict['related_accession'] = record_dict['related_accession'].split(".")[0]
                    if cols[0] == 'gene':
                        std_id = record_dict['chromosome'] + '_' + record_dict['start'] + '_' + record_dict['end'] + '_' + record_dict['strand']
                        gene_std_id[std_id] = record_dict['geneid']
                        gene_dict[record_dict['geneid']] = record_dict
                    elif cols[0] == 'CDS':
                        if 'child_cds' in gene_dict[record_dict['geneid']]:
                            gene_dict[record_dict['geneid']]['child_cds'].append(record_dict)
                        else:
                            gene_dict[record_dict['geneid']]['child_cds'] = [record_dict]

                    else:
                        if 'child' in gene_dict[record_dict['geneid']]:
                            gene_dict[record_dict['geneid']]['child'].append(record_dict)
                        else:
                            gene_dict[record_dict['geneid']]['child'] = [record_dict]
        '''



    def refseq_uniq(self, refseq_dict):
        ## 将基因的下一级为蛋白的加入到基因的下一级
        for gene in refseq_dict:
            if 'child' in refseq_dict[gene]:
                pass
            else:
                refseq_dict[gene]['child'] = []
            children = refseq_dict[gene]['child']
            cds_accs = [child['related_accession'] for child in children]
            if 'child_cds' in refseq_dict[gene]:
                for child_cds in refseq_dict[gene]['child_cds']:
                    if child_cds['product_accession'] in cds_accs:
                        pass
                    else:
                        child_cds['related_accession'] = child_cds['product_accession']
                        child_cds['product_accession'] = ''

                        children.append(child_cds)
        return refseq_dict


    '''
    # 过滤染色体放于各个子数据库里
    def get_genome_info(self, genome_acc):
        print "genome acc is {}".format(genome_acc)
        chr_records = self.genome_api.get_genome_chr(acc_id=genome_acc, assembly_unit={"$in": ["Primary Assembly", "non-nuclear"]})
        chr_list = list()
        for record in chr_records:
            print record
            if record['sequence_role'] == 'assembled-molecule':
                chr_list.append(record['sequence_name'])
            elif record['sequence_role'] in ['unlocalized-scaffold', 'unplaced-scaffold']:
                chr_list.append(record['genbank_accn'])
        return chr_list
    '''

    def filter_ensemble_gene(self, refseq_dict, gene_id2std, ensemble_gene_file, merge_gene_file):
        '''
        refseq_dict : refseq 数据库信息
        gene_id2std :  需要匹配的gene_id 和 std标准对应关系
        '''
        std2gene_dict = dict()
        with open("tmp.log", 'w') as fo:
            fo.write("{}".format(gene_id2std))
        with open(ensemble_gene_file, 'r') as f, open(merge_gene_file, 'w') as fo:
            enterz_list = list()
            std_list = list()
            header = f.readline()
            fields = header.strip("\n").split("\t")
            fo.write("\t".join(fields) + "\n")

            # 写ensemble 已知的基因
            for line in f:
                cols = line.strip("\n").split("\t")
                record_dict = dict(zip(fields, cols))
                sg_gene_id = record_dict['sg_gene_id']
                if sg_gene_id in refseq_dict:
                    record_dict = self.update_byrefseq(record_dict, refseq_dict[sg_gene_id])
                    self.stdgene_dict[sg_gene_id] = record_dict
                    std_list.append(sg_gene_id)
                    self.refgene2stdgene[sg_gene_id] = sg_gene_id
                else:
                    # 只要找到一个id 对应就可以 ？
                    for gid, g2id in self.ref2ens_keys.items():
                        if record_dict[gid] != "\\N":
                            if record_dict[gid] in gene_id2std:
                                std2_id = gene_id2std[record_dict[gid]]
                                self.refgene2stdgene[std2_id] = sg_gene_id
                                record_dict = self.update_byrefseq(record_dict, refseq_dict[std2_id])
                                self.stdgene_dict[sg_gene_id] = record_dict
                                std_list.append(std2_id)
                fo.write("\t".join([record_dict[x] for x in fields]) + "\n")

            # 写不在ensemble 在refseq的基因
            for std2_id in refseq_dict:
                if std2_id in std_list:
                    pass
                else:
                    self.refgene2stdgene[std2_id] = std2_id

                    record_dict = {x:"\\N" for x in fields}
                    record_dict['sg_gene_id'] = std2_id
                    record_dict = self.update_byrefseq(record_dict, refseq_dict[std2_id])
                    self.stdgene_dict[std2_id] = record_dict
                    fo.write("\t".join([record_dict[x] for x in fields]) + "\n")


    def merge_ensemble_refseq_trans(self, tran_dict, refseq_tran_dict, fields):
        '''
        '''
        for gid, t_dic  in tran_dict.items():
            if gid in refseq_tran_dict:
                gene_ref_tran = list()
                ref_tran_dict = refseq_tran_dict[gid]

                # print ref_tran_dict
                for t_id, dic in t_dic.items():
                    if t_id in ref_tran_dict:
                        # 坐标id
                        dic = self.updatetrans_by_refseq(dic, ref_tran_dict[t_id])
                        gene_ref_tran.append(t_id)
                    else:
                        for ref_dict in ref_tran_dict.values():
                            for gid, g2id in self.ref2ens_t_keys.items():
                                if gid in dic and g2id in ref_dict and dic[gid] != "\\N" and dic[gid] == ref_dict[g2id]:
                                    dic = self.updatetrans_by_refseq(dic, ref_dict)
                                    gene_ref_tran.append(t_id)
                                else:
                                    pass
                for ref_t_id, ref_dict in ref_tran_dict.items():
                    if ref_t_id in gene_ref_tran:
                        pass
                    else:
                        dic = {x:"\\N" for x in fields}
                        dic['sg_gene_id'] = gid
                        dic['sg_tran_id'] = ref_t_id
                        # print ref_dict
                        dic = self.updatetrans_by_refseq(dic, ref_dict)
                        ref_tran_dict[ref_t_id] = dic
            else:
                pass
        for ref_g_id in refseq_tran_dict:
            if ref_g_id not in tran_dict:
                tran_dict[ref_g_id] = dict()
                for ref_t_id,ref_dic in refseq_tran_dict[ref_g_id].items():
                    dic = {x:"\\N" for x in fields}
                    dic['sg_gene_id'] = ref_g_id
                    dic['sg_tran_id'] = ref_t_id
                    dic = self.updatetrans_by_refseq(dic, ref_dict)
                    tran_dict[ref_g_id][ref_t_id] = dic


    def filter_ensemble_tran(self, refseq_tran_dict, tran_id2std, ensemble_tran_file, merge_tran_file):
        '''
        refseq_dict : refseq 数据库信息
        tran_id2std :  需要匹配的tran_id 和 std标准对应关系
        '''
        std2tran_dict = dict()

        tran_dict = dict()

        with open(ensemble_tran_file, 'r') as f:
            header = f.readline()
            fields = header.strip("\n").split("\t")

        with open(ensemble_tran_file, 'r') as f:
            for dic in csv.DictReader(f, delimiter='\t'):
                if dic["sg_gene_id"] in tran_dict:
                    tran_dict[dic["sg_gene_id"]].update({dic["sg_tran_id"]: dic})
                else:
                    tran_dict[dic["sg_gene_id"]] ={dic["sg_tran_id"]: dic}

        self.merge_ensemble_refseq_trans(tran_dict, refseq_tran_dict, fields)

        with open(merge_tran_file, 'w') as fo:
            fo.write("\t".join(fields) + "\n")
            for g, t_dic in tran_dict.items():
                for t, record_dict in t_dic.items():
                    fo.write("\t".join([record_dict[x] for x in fields]) + "\n")


    def merge_ensemble_gene_trans(self, ens_gid, trans_list, fo, fields):

        print "merge gene {}".format(ens_gid)
        if ens_gid in self.gene2std:
            std_id = self.gene2std[ens_gid]
        else:
            # 基因在其它染色体上被过滤掉
            return
        # gene_list = self.std2gene[std_id]
        ref_trans_list = list()
        if std_id in self.std2enterz:
            ref_trans_list = self.refseq_dict[self.std2enterz[std_id]]['child']
        if len(ref_trans_list) == 0:
            pass
        else:
            matched_ref_trans = list()
            for trans in trans_list:
                for n, a_trans in enumerate(ref_trans_list):
                    ##
                    match1 = [trans[x] for x in self.trans_key_matching.keys()]
                    match2 = [a_trans[x] for x in self.trans_key_matching.values()]
                    if len(set(match1).intersection(set(match2))) > 0:
                        trans = self.updatetrans_by_refseq(trans, a_trans)
                        matched_ref_trans.append(n)
                        fo.write("\t".join([std_id] + [trans[x] for x in fields]) + "\n")
                        break
            # 合并没有匹配到的refseq转录本
            for n, a_trans in enumerate(ref_trans_list):
                if n in matched_ref_trans:
                    pass
                else:
                    trans = {x:"\\N" for x in fields}
                    trans = self.updatetrans_by_refseq(trans, a_trans)
                    fo.write("\t".join([std_id] + [trans[x] for x in fields]) + "\n")


    '''
    def filter_ensemble_trans(self, refseq_dict, std_gene_id, ensemble_tran_file, merge_tran_file):
        self.refseq_dict = self.refseq_uniq(refseq_dict)

        with open(ensemble_tran_file, 'r') as f, open(merge_tran_file, 'w') as fo:
            enterz_list = list()
            header = f.readline()
            fields = header.strip("\n").split("\t")
            fo.write("\t".join(["sg_gene_id"] + fields) + "\n")

            origin_gene_id = ""
            # 写ensemble 已知的转录本
            tran_list = []

            for line in f:
                cols = line.strip("\n").split("\t")
                record_dict = dict(zip(fields, cols))

                gene_id = record_dict['gene_id_1020_key']
                if gene_id != origin_gene_id and origin_gene_id != "":
                    self.merge_ensemble_gene_trans(origin_gene_id, tran_list, fo, fields)
                    trans_list = []

                origin_gene_id = gene_id
                tran_list.append(record_dict)

            # 最后一条记录
            self.merge_ensemble_gene_trans(origin_gene_id, tran_list, fo, fields)
            for std_id,eid in std_gene_id.items():
                if std_id in self.use_enterz_genes:
                    ref_trans_list = self.refseq_dict[eid]['child']
                    for a_trans in ref_trans_list:
                        trans = {x:"\\N" for x in fields}
                        trans = self.updatetrans_by_refseq(trans, a_trans)
                        fo.write("\t".join([std_id] + [trans[x] for x in fields]) + "\n")

    '''
    def update_byrefseq(self, origin_dict, refseq_dict):
        for field, ref_field in self.use_refseq_fields.items():
            if refseq_dict[ref_field] != "\\N" and refseq_dict[ref_field] != "":
                origin_dict[field] = refseq_dict[ref_field]
        return origin_dict

    def updatetrans_by_refseq(self, origin_dict, refseq_dict):
        for field, ref_field in self.use_refseq_t_fields.items():
            if refseq_dict[ref_field] != "\\N" and refseq_dict[ref_field] != "":
                origin_dict[field] = refseq_dict[ref_field]
        return origin_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-acc', type=str, required=False,
                        help='genome acc')
    parser.add_argument('-gene_table', type=str, required=True,
                        help='refseq gene table')
    parser.add_argument('-tran_table', type=str, required=True,
                        help='refseq tran table')
    parser.add_argument('-ensemble_gene', type=str, default=True,
                        help="ensemble gene merged tsv file")
    parser.add_argument('-ensemble_tran', type=str, default=True,
                        help="ensemble tran merged tsv file")
    # ----------------------------------------------------------------------------------------------
    args = parser.parse_args()
    mg = MergeGene(None)
    mg.merge_pipe(args.gene_table, args.tran_table, args.ensemble_gene, args.ensemble_tran)


