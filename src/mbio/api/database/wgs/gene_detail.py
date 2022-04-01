# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.05.07

from api_base import ApiBase
from bson.objectid import ObjectId
import datetime
import os
import re


class GeneDetail(ApiBase):
    def __init__(self, bind_object):
        """
        基因详情页导表
        """
        super(GeneDetail, self).__init__(bind_object)
        self._project_type = "dna_wgs"

    def add_sg_gene_sample_seq_detail(self, seq_id, gene_id, gene_path):
        """
        sg_sample_seq_detail
        gene_path: gene_id 对应的样本序列表
        """
        seq_id = self.check_objectid(seq_id)
        self.check_exists(gene_path)
        data_list = []
        with open(gene_path, "r") as f:
            lines = f.readlines()
            sum = len(lines)
            for i in range(sum-1):
                m = re.match(r">(.+)", lines[i])
                if m:
                    specimen_id = m.group(1)
                    seq = lines[i+1].strip()
                    insert_data = {
                        "gene_sample_id": seq_id,
                        "gene_id": gene_id,
                        "specimen_id": specimen_id,
                        "length": len(seq),
                        "seq": seq
                    }
                    data_list.append(insert_data)
        if not data_list:
            self.bind_object.logger.info("没有找到对应的样本序列")
            # raise Exception("没有找到对应的样本序列")
        else:
            self.col_insert_data("sg_gene_seq_detail", data_list)

    def add_sg_gene_structure_detail(self, structure_id, gene_id, start, end, exon_path, alt_path, mirna_path):
        """
        sg_gene_structure_detail
        基因结构导表
        exon_path: 外显子序列及种类
        alt_path：变异位点信息
        mirna_path：mirna的其实终止位置
        >chr1   rna1577 id10216 8281766 8281896
        """
        structure_id = self.check_objectid(structure_id)
        self.check_exists(exon_path)
        self.check_exists(alt_path)
        self.check_exists(mirna_path)
        seq_data_list, graph_data_list, categories = [], [], []
        graph_dict, seq_dict, loca_dict = {}, {}, {}
        with open(exon_path, "r") as f:
            lines = f.readlines()
            for i in range(len(lines)-1):
                if lines[i].startswith(">"):                                # startwith">"行处理！
                    # item = lines[i].strip().split("\t")
                    # type = item[1]                                          # 第二列存了type！！！！
                    # if type not in categories:
                    #     categories.append(type)                             # type存入categories.append(type)
                    #     graph_dict[type] = {}                               # graph_dict[type] = {}
                    #     seq_dict[type] = lines[i+1].strip()                 # seq_dict[type] = lines[i+1].strip() 存入id的seq
                    # graph_dict[type][item[2]] = [int(item[3]), int(item[4])]    # 存入rna1577 (id10216) = start end  8281766 8281896
                    item = lines[i].strip(">")
                    item = item.strip()
                    seq_dict[item] = lines[i+1].strip()     
        with open(mirna_path, "r") as f:
            lines = f.readlines()
            for line in lines:
                item = line.strip().split("\t")
                if item[2] not in loca_dict.keys():
                    loca_dict[item[2]]=[]                
                if (item[2] not in graph_dict.keys()):      # 只有type（rna1577）
                    graph_dict[item[2]] = {} 
                graph_dict[item[2]][item[3]] = item[4]
                loca_dict[item[2]].append(int(item[3]))
                loca_dict[item[2]].append(int(item[4]))                
                # loca_dict[item[2]] = [int(item[3]), int(item[4])]           #item[2]是type:{rna0};       eg: >chr1  exon    {rna0}    {4706    5095} # 覆盖
                if item[2] not in categories:
                    categories.append(item[2])
        snp_var_dict, indel_var_dict = {}, {}
        with open(alt_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                v = item[2] + "," + item[3]                                 # #CHROM    POS REF ALT
                var_type = True
                for k in v.split(","):
                    if len(k) != 1:
                        var_type = False
                if var_type:
                    snp_var_dict[item[1]] = [int(item[1]), item[2], item[3]]
                else:
                    indel_var_dict[item[1]] = [int(item[1]), item[2], item[3]]
        snp_var_pos = snp_var_dict.keys()
        indel_var_pos = indel_var_dict.keys()
        insert_data = {
            "origin_id": structure_id,
            "name": "基因结构图",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "categories": categories,
            "legend": ["exon"],
            "location": "gene_structure",
            "type": 1,
            "x_start": start,
            "x_end": end
        }
        main_id = self.db["sg_structure"].insert_one(insert_data).inserted_id
        for name in graph_dict.keys():
            types, values, snp_datas, indel_datas = [], [], [], []
            for i in graph_dict[name].keys():
                types.append("exon")
            for start_pos in graph_dict[name]:
                a = [int(start_pos), int(graph_dict[name][start_pos])]
                values.append(a)
        # for name in categories:
        #     types = graph_dict[name].keys()
        #     type, values, snp_datas, indel_datas = [], [], [], []
        #     for t1 in types:
        #         type.append("exon")
        #         values.append(graph_dict[name][t1])
            # for pos in snp_var_pos:
            #     if loca_dict[name][0] <= int(pos) <= loca_dict[name][1]:
            #         snp_datas.append(snp_var_dict[pos])
            # for pos in indel_var_pos:
            #     if loca_dict[name][0] <= int(pos) <= loca_dict[name][1]:
            #         indel_datas.append(indel_var_dict[pos])
            for pos in snp_var_pos:
                if min(loca_dict[name]) <= int(pos) <= max(loca_dict[name]):
                    snp_datas.append(snp_var_dict[pos])
            for pos in indel_var_pos:
                if min(loca_dict[name]) <= int(pos) <= max(loca_dict[name]):
                    indel_datas.append(indel_var_dict[pos])    
            seq_data = {
                "structure_id": structure_id,
                "gene_id": gene_id,
                "name": name,
                "len": len(seq_dict[name]),
                # "location": gene_id + ":" + str(loca_dict[name][0]) + "-" + str(loca_dict[name][1]),
                "location": gene_id + ":" + bytes(start) + "-" + bytes(end),
                "seq": seq_dict[name]
            }
            seq_data_list.append(seq_data)
            graph_data = {
                "structure_id": main_id,
                "name": name,
                "type": types,
                "values": values,
                "snp_datas": snp_datas,
                "indel_datas": indel_datas
            }
            graph_data_list.append(graph_data)
        if seq_data_list:
            self.col_insert_data("sg_gene_structure_detail", seq_data_list)
            self.col_insert_data("sg_structure_detail", graph_data_list)
            self.update_db_record("sg_gene_structure", {"_id": structure_id}, {"categories": categories})
        else:
            self.bind_object.logger.info("在此基因{}上没找到转录本".format(gene_id))


if __name__ == "__main__":
    a = GeneDetail(None)
    seq_id = "5aefc82aa4e1af6fae4e4cf3"
    gene_id = "gene0"
    gene_path = "/mnt/ilustre/users/sanger-dev/workspace/20180507/Single_gene_samples_seq/GeneSamplesSeq/output/gene0.fa"
    # a.add_sg_gene_sample_seq_detail(seq_id, gene_id, gene_path)
    structure_id = "5af267ada4e1af50802eb333"
    exon_path = "/mnt/ilustre/users/sanger-dev/sg-users/cuiqingmei/Single_tsanger_30192_0608085210694902_2021/GeneStructure/output/gene10315.exon.fa"
    alt_path = "/mnt/ilustre/users/sanger-dev/sg-users/cuiqingmei/Single_tsanger_30192_0608085210694902_2021/GeneStructure/output/gene10315.xls"
    mirna_path = "/mnt/ilustre/users/sanger-dev/sg-users/cuiqingmei/Single_tsanger_30192_0608085210694902_2021/GeneStructure/output/gene10315.mirna.xls"
    a.add_sg_gene_structure_detail(structure_id, gene_id, 1, 10000, exon_path, alt_path, mirna_path)
