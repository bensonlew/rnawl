# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# modified 2018.0825

from api_base import ApiBase
from collections import defaultdict
import os
import datetime
import json


class GeneSampleSeq(ApiBase):
    """
    群体进化，GWAS关联分析接口
    """
    def __init__(self, bind_object):
        super(GeneSampleSeq, self).__init__(bind_object)
        self._project_type = "dna_evolution"

    def add_sg_structure(self, gene_id, structure_id, exon_path):
        """
        基因结构图导表
        sg_structure
        sg_structure_detail
        structure_id:上一级主表sg_gene_structure id
        """
        self.check_exists(exon_path)
        start_list = []
        end_list = []
        legend_list = []
        categories = []
        graph_data_list = []
        all_data = []
        insert_data_list = []
        with open(exon_path, "r") as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith(">"):
                    line_strip = line.strip(">")
                    line_strip_blank = line_strip.strip()
                    line_split = line_strip_blank.split("\t")
                    start_list.append(line_split[-3])
                    end_list.append(line_split[-2])
                    legend_list.append(line_split[-1])
                    categories.append(line_split[-5])
                    detail_dict = {}
                    detail_dict["name"] = line_split[-5]
                    detail_dict["start"] = line_split[-3]
                    detail_dict["end"] = line_split[-2]
                    detail_dict["type"] = line_split[-1]
                    all_data.append(detail_dict)
        start = min(start_list)
        end = max(end_list)
        legend_list = list(set(legend_list))
        categories = list(set(categories))
        origin_id = self.check_objectid(structure_id)
        insert_data = {
            "origin_id": origin_id,
            "name": "基因结构图",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "categories": categories,
            "legend": legend_list,
            "location": "gene_structure",
            "type": 1,
            "x_start": int(start),
            "x_end": int(end)
        }
        insert_data_list.append(insert_data)
        structure_id = self.col_insert_data("sg_structure", insert_data_list)
        for name in categories:
            values = []
            types = []
            for data_dict in all_data:
                if data_dict["name"] == name and data_dict["type"] == "exon":
                    data_list = []
                    data_list.append(int(data_dict["start"]))
                    data_list.append(int(data_dict["end"]))
                    values.append(data_list)
                    types.append(data_dict["type"])
            for data_dict in all_data:
                if data_dict["name"] == name and data_dict["type"] == "CDS":
                    data_list = []
                    data_list.append(int(data_dict["start"]))
                    data_list.append(int(data_dict["end"]))
                    values.append(data_list)
                    types.append(data_dict["type"])
            graph_data = {
                "structure_id": structure_id,
                "name": name,
                "type": types,
                "values": values,
                "snp_datas": [],
                "indel_datas": []
            }
            graph_data_list.append(graph_data)
        if graph_data_list:
            self.col_insert_data("sg_structure_detail", graph_data_list)
            # self.update_db_record("sg_gene_structure", {"_id": structure_id}, {"categories": categories})
        else:
            self.bind_object.logger.info("在此基因{}上没找到转录本".format(gene_id))

    def add_sg_gene_sample_seq_detail(self, seq_id, gene_id, gene_path, type, status, desc):
        """
        sg_sample_seq_detail
        gene_path: gene_id 对应的样本序列表
        """
        seq_id = self.check_objectid(seq_id)
        # self.check_exists(gene_path)
        insert_data = {
            "origin_id": seq_id,
            "gene_id": gene_id,
            "seq": gene_path,
            "type": type,
            "status": status,
            "desc": desc
        }
        self.db["sg_gene_seq_detail"].insert_one(insert_data)

    def check_exists(self, file_path):
        """
        用于检查文件及文件夹是否存在
        :param file_path:
        :return:
        """
        if not os.path.exists(file_path):
            raise Exception("文件或文件夹{}不存在！".format(file_path))


if __name__ == "__main__":
    a = GeneSampleSeq(None)
    # a.add_sg_gene_sample_seq_detail("5b6fb098f6b9e414d0f4bca0", "genemark-tig00000001-processed-gene-0.11",
    #                                 "/mnt/ilustre/users/sanger-dev/workspace/20180927/Single_test_gene_sample_seq_0927_1/GeneSampleSeq/output/protein.fa", "protein")
    a.add_sg_structure("gene0", "5b6fb098f6b9e414d0f4bca0", "/mnt/ilustre/users/sanger-dev/workspace/20180927/Single_test_gene_sample_seq_0927_1/GeneSampleSeq/output/gene0.exon.fa")
