# -*- coding: utf-8 -*-
import os
import math
import numpy as np
import argparse
import glob
import re
import gzip
import json
import copy
import csv
from collections import defaultdict
from collections import OrderedDict
# from mbio.api.database.gene_db import genome
from mbio.api.database.gene_db import logic

class MergeUniprot(object):
    def __init__(self):
        # self.genome_api = genome.Genome(None)
        # self.genome_acc = genome_acc
        self.logic_api = logic.Logic(None)
        self.uniprot_dict = dict()
        self.uniprot_tran_dict = dict()
        # self.ref2ens_keys = OrderedDict({})

        self.uniprot_match_tran = OrderedDict({
            "Ensembl_TRS": "ensembl_transcript_id",
            "Ensembl_PRO": "ensembl_peptide_id",
            "RefSeq": "refseq_peptide",
            "RefSeq_NT": "refseq_ncrna"
        })

        # 后续需要合并refseq_mrna refseq_dna refseq_ncrna

        self.uniprot_match = self.uniprot_match_tran

        # 由其它ID对应到uniprot ID
        self.uniprot_match_ids = OrderedDict({})
        self.uniprot_match_tran_ids = OrderedDict({})

        for key in self.uniprot_match.keys():
            self.uniprot_match_ids[key] = dict()
            self.uniprot_match_tran_ids[key] = dict()


    def logic_api_info(self):
        record_list = self.logic_api.get_logic_base()
        self.use_uniprot_dict = dict()
        self.has_uniprot_list = list()
        self.append_uniprot_list = list()

        for record in record_list:
            if record.get("UNIPROT", "") != "":
                self.has_uniprot_list.append(record.get("UNIPROT"))
            if record.get("UNIPROT", "") != "" and record.get("uniprot_use", 0.0) != 0.0:
                self.use_uniprot_dict[record.get("UNIPROT", "")] = record.get("field", "")

        self.add_fields = self.use_uniprot_dict.values()
        self.matching_list = list()

        '''
        for record in record_list:
            if record.get("uniprot_match", "") != "":
                self.matching_list.append((record["uniprot_match"], [record["field"], record["UNIPROT"]]))

        self.matching_dict = OrderedDict(sorted(self.matching_list, key=lambda x: x[1][0]))
        self.matching_dict
        '''
        self.matching_dict = OrderedDict(sorted(self.matching_list, key=lambda x: x[1][0]))

    def read_idmapping_list(self, id_file):
        with gzip.open(id_file, "r") as f:
            for line in f:
                uni_id, db, db_id = line.strip().split("\t")
                if db == "EnsemblGenome":
                    db = "Ensembl"
                if db == "EnsemblGenome_TRS":
                    db = "Ensembl_TRS"
                if db == "EnsemblGenome_PRO":
                    db = "Ensembl_PRO"
                if db not in self.has_uniprot_list and db not in self.append_uniprot_list:
                    self.append_uniprot_list.append(db)
                if len(uni_id) == 6 or uni_id[6] != "-":
                    if uni_id in self.uniprot_dict:
                        if db in self.uniprot_dict[uni_id]:
                            self.uniprot_dict[uni_id][db].add(db_id)
                        else:
                            self.uniprot_dict[uni_id][db] = set([db_id])
                    else:
                        self.uniprot_dict[uni_id] = {db: set([db_id])}

                    if db in self.uniprot_match:
                        ## 出现refseq对应多个uniprot的情况
                        ## Q9M0T4  RefSeq_NT       NM_202792.1
                        ## A0A384KHN3      RefSeq_NT       NM_202792.1
                        if db_id.split(".")[0] in self.uniprot_match_ids[db]:
                            pass
                        else:
                            self.uniprot_match_ids[db][db_id.split(".")[0]] = uni_id

                elif uni_id[6] == "-":
                    # 转录本水平的映射
                    if uni_id in self.uniprot_tran_dict:
                        if db in self.uniprot_tran_dict[uni_id]:
                            self.uniprot_tran_dict[uni_id][db].add(db_id)
                        else:
                            self.uniprot_tran_dict[uni_id][db] = set([db_id])
                    else:
                        self.uniprot_tran_dict[uni_id] = {db: set([db_id])}

                    if db in self.uniprot_match_tran:
                        self.uniprot_match_tran_ids[db][db_id.split(".")[0]] = uni_id



        # uniprot 数据添加额外的字段
        for append_id in self.append_uniprot_list:
            self.use_uniprot_dict.update({append_id: append_id.lower()})
        self.add_fields = self.use_uniprot_dict.values()

        # print self.uniprot_dict

    def merge_pipe(self, id_list, merge_tran, merge_tran_out):
        self.logic_api_info()
        self.read_idmapping_list(id_list)
        self.merge_uniprot_trans(merge_tran, merge_tran_out)

    def merge_gene_dic(self, iso_id, dic):
        '''
        iso_id: 转录本id
        dic: 转录本属性
        '''
        # 合并蛋白信息到转录本上

        # print "dic is {}".format(dic)
        # print "uni dic is {}".format( self.uniprot_dict )
        uni_id = iso_id.split("-")[0]
        if uni_id  in self.uniprot_dict:
            uniprot_dict1 = self.uniprot_dict[uni_id]
            # 加[] 会不会更好
            dic['UniProtKB-ID'] = uni_id
            for key in uniprot_dict1:
                if key in dic:
                    pass
                else:
                    dic[key] = uniprot_dict1[key]
        else:
            pass
        return dic

    def merge_uniprot_dict(self, dic, uniprot_dic1):
        # 合并uniprot 到 原始结果
        # print "uniprot dict is {}".format(uniprot_dic1)
        for u_key, m_key in self.use_uniprot_dict.items():
            # print u_key, m_key
            if m_key in dic:
                if dic[m_key] in ["\\N", ""]:
                    if type(uniprot_dic1.get(u_key, "")) == str:
                        dic[m_key] = uniprot_dic1.get(u_key, "")
                    else:
                        dic[m_key] = "|".join(uniprot_dic1.get(u_key, ""))
                else:
                    pass
            else:
                if type(uniprot_dic1.get(u_key, "")) == str:
                    dic[m_key] = uniprot_dic1.get(u_key, "")
                else:
                    dic[m_key] = "|".join(uniprot_dic1.get(u_key, ""))
        return dic

    def merge(self, dic):
        # 优先映射转录本

        for u_key, m_key in self.uniprot_match_tran.items():
            if m_key in dic and dic[m_key] != "" and dic[m_key] in self.uniprot_match_tran_ids[u_key]:
                dic['uniprot_isoform'] = self.uniprot_match_tran_ids[u_key][dic[m_key]]
                iso_id = dic['uniprot_isoform']
                uniprot_dic_1 = self.merge_gene_dic(iso_id, self.uniprot_tran_dict[dic['uniprot_isoform']])
                # print uniprot_dic_1
                # print dic
                dic = self.merge_uniprot_dict(dic, uniprot_dic_1)
                # print dic
                return dic

        # 根据UNIPROT ID 映射
        # 转录本，基因用 uniprot_gn_id
        if dic['uniprotswissprot'] in self.uniprot_dict and dic['uniprotswissprot'] != "\\N":
            dic['uniprot_gn_id'] = dic['uniprotswissprot']
            dic = self.merge_uniprot_dict(dic, self.uniprot_dict[dic['uniprotswissprot']])
            return dic
        elif dic.get('uniprot_gn_id', "\\N") in self.uniprot_dict:
            dic = self.merge_uniprot_dict(dic, self.uniprot_dict[dic['uniprot_gn_id']])
            return dic
        elif dic.get('uniprotsptrembl', "\\N") in self.uniprot_dict:
            dic['uniprot_gn_id'] = dic['uniprotsptrembl']
            dic = self.merge_uniprot_dict(dic, self.uniprot_dict[dic['uniprot_gn_id']])
            return dic

        # 更具其它id映射
        for u_key, m_key in self.uniprot_match.items():
            if m_key in dic and dic[m_key] != "" and dic[m_key] in self.uniprot_match_ids[u_key]:
                dic['UniProtKB-ID'] = self.uniprot_match_ids[u_key][dic[m_key]]
                uniprot_dic_1 = self.uniprot_dict[dic['UniProtKB-ID']]
                dic = self.merge_uniprot_dict(dic, uniprot_dic_1)
                return dic

        return dic

    def merge_uniprot_trans(self, merge_tran, merge_tran_out):
        '''
        merge_tran : refseq 数据库信息
        merge_tran_out :  需要匹配的tran_id 和 std标准对应关系
        '''

        with open(merge_tran, 'r') as f:
            header = f.readline()
            fields = header.strip("\n").split("\t")

        if not ("uniprot_gn_id" in fields):
            fields.append("uniprot_gn_id")
        fields += ["uniprot_isoform"]
        for field in self.add_fields:
            if field in fields:
                pass
            else:
                fields.append(field)
        with open(merge_tran, 'r') as f, open(merge_tran_out, 'w') as fo:
            fo.write("\t".join(fields) + "\n")
            for dic in csv.DictReader(f, delimiter='\t'):
                dic = self.merge(dic)
                fo.write("\t".join([dic.get(x, "\\N") for x in fields]) + "\n")
        with open(merge_tran_out + ".fields", 'w') as fo:
            for k,v in self.use_uniprot_dict.items():
                fo.write("{}\t{}\n".format(k, v))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-uniprot', type=str, required=True,
                        help='uniprot list ')
    parser.add_argument('-merge_in', type=str, required=True,
                        help='input merged transcript file')
    parser.add_argument('-merge_out', type=str, default=True,
                        help="output merged transcript file")

    args = parser.parse_args()
    mg = MergeUniprot()
    mg.merge_pipe(args.uniprot, args.merge_in, args.merge_out)


