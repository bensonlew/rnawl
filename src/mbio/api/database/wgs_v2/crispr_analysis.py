# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 20190304

from api_base import ApiBase
import datetime
from collections import defaultdict
import re


class CrisprAnalysis(ApiBase):
    def __init__(self, bind_object):
        """
        crispr_analysis接口
        """
        super(CrisprAnalysis, self).__init__(bind_object)
        self._project_type = "dna_wgs_v2"

    def sg_crispr_analysis(self, task_id, project_sn, params):
        """
        变异位点比较分析主表
        :return:
        """
        main_id = self.add_main_table("sg_crispr_analysis", task_id, project_sn, params, "crispr_analysis",
                                      "基因编辑接口导表", "开始进行CRISPR分析导表")
        self.update_db_record("sg_crispr_analysis", {"_id": main_id}, {"main_id": main_id, "status": "end"})
        return main_id

    def add_modified_effiency(self, file_path, gene_name, bam_list, crispr_id, task_id):
        """"
        用户传入的位置内可能有多处被编辑，但是每个样本没处的可能性只有两个，也就是crispr只应用于二倍体
        """
        data_list = []
        crispr_id = self.check_objectid(crispr_id)
        with open(file_path)as f:
            pos_list = []  # 用户输入的位置并不是一定只有一个敲除位置，也就是说不一定对应的只有一个R。后来确定的只有一种情况，但是为了修改方便，依然保留list。
            lines = f.readlines()
            for line in lines[1:]:
                pos = line.strip().split("\t")[0]
                if pos not in pos_list:
                    pos_list.append(pos)
        with open(file_path) as n:
            lines = n.readlines()
            for item in pos_list:
                modified_num = 0
                variant_dict = {}
                for line in lines[1:]:
                    if line.strip().split("\t")[0] == item:
                        if line.strip().split("\t")[1] not in variant_dict.keys():
                            variant_dict[line.strip().split("\t")[1]] = []
                        if len(variant_dict[line.strip().split("\t")[1]]) < 3:
                            variant_dict[line.strip().split("\t")[1]].append(line.strip().split("\t")[4])  # 这个字典的key为样本名，value为基因型
                    else:
                        continue
                print str(variant_dict) + "&&&&&&&&&&&&&&&&&&&&&&&&&&&^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
                for sample in bam_list:
                    if len(variant_dict[sample]) == 1:
                        if variant_dict[sample].count("wt") == 0:
                            modified_num += 2
                    elif len(variant_dict[sample]) == 2:
                        if variant_dict[sample].count("wt") == 0:
                            modified_num += 2
                        elif variant_dict[sample].count("wt") == 1:
                            modified_num += 1
                modified_effiency = str(float(modified_num) * 100 / (2 * len(bam_list))) + "%"
                gene_name1 = gene_name + " " + item if gene_name != "" else item
                insert_data = {
                    "target_gene": gene_name1,
                    "modified": modified_num,
                    "modified_effiency": modified_effiency,
                    "crispr_id": crispr_id,
                    "task_id": task_id
                }
                data_list.append(insert_data)
        if len(data_list) == 0:
            self.bind_object.logger.info("{}文件为空！".format(file_path))
        else:
            self.col_insert_data("sg_crispr_stat", data_list)

    def add_crispr_stat(self, file_path, crispr_id, task_id):
        crispr_id = self.check_objectid(crispr_id)
        with open(file_path) as f:
            sum1 = 0  # 用于统计allele的总数
            i = 0  # 用于统计插入的个数
            d = 0
            c = 0
            s = 0
            wt = {}  # 因为wt需要写序列，所以wt以这种形式字典形式来存储。
            data_list = []
            sequence_list = []
            mutation_type_list = []
            frequence_list = []
            lines = f.readlines()
            """
            后续修改各个样本也要存为list的形式，因为wt为之前所写，因此保持不变"""
            sequence_dict = {}
            mutation_type_dict = {}
            allele_freq_dict = {}
            for line in lines[1:]:
                sum1 += 1
                if line.strip().split("\t")[4] == "i":
                    i += 1
                elif line.strip().split("\t")[4] == "d":
                    d += 1
                elif line.strip().split("\t")[4] == "c":
                    c += 1
                elif line.strip().split("\t")[4] == "s":
                    s += 1
                elif line.strip().split("\t")[4] == "wt":
                    if line.strip().split("\t")[3] not in wt.keys():
                        wt[line.strip().split("\t")[3]] = 1
                    else:
                        wt[line.strip().split("\t")[3]] += 1
                if re.match("new.*", line.strip().split("\t")[1]):
                    pass
                else:
                    # insert_data = {
                    #     "sample_id":  line.strip().split("\t")[1],
                    #     "target_seq": line.strip().split("\t")[3],
                    #     "mutation_type": line.strip().split("\t")[4] + line.strip().split("\t")[5]
                    #     if line.strip().split("\t")[4] == "i" or line.strip().split("\t")[4] == "d" else line.strip().split("\t")[4],
                    #     "allele_freq": "-",
                    #     "crispr_id": crispr_id,
                    #     "task_id": task_id
                    # }
                    # data_list.append(insert_data)
                    if line.strip().split("\t")[1] not in sequence_dict.keys():
                        sequence_dict[line.strip().split("\t")[1]] = []
                    sequence_dict[line.strip().split("\t")[1]].append(line.strip().split("\t")[3])
                    if line.strip().split("\t")[1] not in mutation_type_dict.keys():
                        mutation_type_dict[line.strip().split("\t")[1]] = []
                        allele_freq_dict[line.strip().split("\t")[1]] = []
                    mutation_type_dict[line.strip().split("\t")[1]].append(line.strip().split("\t")[4] + line.strip().
                                                                           split("\t")[5] if
                                                                           line.strip().split("\t")[4] == "i" or line.
                                                                           strip().
                                                                           split("\t")[4] == "d" else line.strip().
                                                                           split("\t")[4])
                    allele_freq_dict[line.strip().split("\t")[1]].append("-")
        for keys in sequence_dict.keys():
            insert_data = {
                "sample_id": keys,
                "target_seq": sequence_dict[keys],
                "mutation_type": mutation_type_dict[keys],
                "allele_freq": allele_freq_dict[keys],
                "task_id": task_id,
                "crispr_id": crispr_id,
            }
            data_list.append(insert_data)
        for keys in wt.keys():
            sequence_list.append(keys)
            mutation_type_list.append("wt")
            frequence_list.append(str(wt[keys]) + "/" + str(sum1))
        for m in range(4):
            sequence_list.append("********************")
        mutation_type_list.append("i")
        frequence_list.append(str(i) + "/" + str(sum1))
        mutation_type_list.append("d")
        frequence_list.append(str(d) + "/" + str(sum1))
        mutation_type_list.append("s")
        frequence_list.append(str(s) + "/" + str(sum1))
        mutation_type_list.append("c")
        frequence_list.append(str(c) + "/" + str(sum1))
        insert_data1 = {
            "sample_id": "All",
            "target_seq": sequence_list,
            "mutation_type": mutation_type_list,
            "allele_freq": frequence_list,
            "crispr_id": crispr_id,
            "task_id": task_id
        }
        data_list.insert(0, insert_data1)
        if len(data_list) == 0:
            self.bind_object.logger.info("{}文件为空！".format(file_path))
        else:
            self.col_insert_data("sg_crispr_stat_detail", data_list)
        x_category = ["i", "d", "s", "c"]
        bar_id = self.add_sg_bar(crispr_id, task_id, x_category, "crispr_type_distribution")
        self.add_sg_bar_detail(bar_id, "crispr_type_distribution", [i, d, s, c])

    def add_crispr_off_target(self, file_path, crispr_id, task_id):
        crispr_id = self.check_objectid(crispr_id)
        with open(file_path) as f:
            lines = f.readlines()
            data_list = []
            n_modified = 0
            total = 0
            for line in lines[1:]:
                n_modified += int(line.strip().split("\t")[1])
                total += int(line.strip().split("\t")[2])
                insert_data = {
                    "sample_name": line.strip().split("\t")[0],
                    "putative_off_target": str(line.strip().split("\t")[1]) + "/" + str(line.strip().split("\t")[2]),
                    "off_target_rate": str(round(float(line.strip().split("\t")[1]) * 100 / int(line.strip().split("\t")
                                                                                                [2]), 2)) + "%",
                    "crispr_id": crispr_id,
                    "task_id": task_id,
                }
                data_list.append(insert_data)
            insert_data = {
                "sample_name": "All",
                "putative_off_target" : str(n_modified) + "/" + str(total),
                "off_target_rate": str(round(float(n_modified)*100/total, 2)) + "%",
                "crispr_id": crispr_id,
                "task_id": task_id,

            }
            data_list.insert(0, insert_data)
            if len(data_list) == 0:
                self.bind_object.logger.info("{}文件为空！".format(file_path))
            else:
                self.col_insert_data("sg_crispr_off_target", data_list)

    def add_sg_bar(self, feature_id, task_id, categories, location):
        """
        sg_bar
        """
        insert_data = {
            "origin_id": feature_id,
            "task_id": task_id,
            "location": location,
            "type": 1,
            "categories": categories,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "other_attr": "",
            "ext_name": ""
        }
        bar_id = self.db["sg_bar"].insert_one(insert_data).inserted_id
        return bar_id

    def add_sg_bar_detail(self, bar_id, name, values):
        """
        sg_bar_detail
        """
        bar_id = self.check_objectid(bar_id)
        insert_data = {
            "bar_id": bar_id,
            "name": name,
            "value": values
        }
        self.db["sg_bar_detail"].insert_one(insert_data)


if __name__ == "__main__":
    a = CrisprAnalysis(None)
    task_id = "wgs_v2_test"
    project_sn = "wgs_v2"
    name = "crispr_test"
    crispr_id = "5cad77be17b2bf11dc05bff2"
    # a.sg_crispr_analysis(task_id, project_sn, "为crispr创造主表")
    # a.add_crispr_off_target("/mnt/ilustre/users/sanger-dev/workspace/20190404/Single_test_crispr_module20190404100108/"
    #                         "CrisprAnalysis/output/off_target.txt", crispr_id, task_id)
    # a.add_modified_effiency("/mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/test/on_target.txt", "p53", ["S2_12",
    #                                                                                                         "S2_5",
    #                                                                                                         "S2_6"],
    #                         crispr_id, task_id)
    a.add_crispr_stat("/mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/test/on_target.txt", crispr_id, task_id)



