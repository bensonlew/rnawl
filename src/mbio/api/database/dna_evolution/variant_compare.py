# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 20180827

from api_base import ApiBase
import datetime
from collections import defaultdict
import re


class VariantCompare(ApiBase):
    def __init__(self, bind_object):
        """
        比较分析接口的所有导表
        """
        super(VariantCompare, self).__init__(bind_object)
        self._project_type = "dna_evolution"

    def sg_variant_compare(self, task_id, project_sn, params):
        """
        变异位点比较分析主表
        :return:
        """
        main_id = self.add_main_table("sg_variant_compare", task_id, project_sn, params, "origin_variant_compare",
                                      "比较分析接口导表", "开始进行变异位点比较插图")
        self.update_db_record("sg_variant_compare", {"_id": main_id}, {"main_id": main_id, "status": "end"})
        return main_id

    def update_varaint_compare(self, download_path, filter_vcf_path, update_id):
        update_id = self.check_objectid(update_id)
        self.update_db_record("sg_variant_compare", {"_id": update_id}, {"pop_table_path": download_path,
                                                                         "vcf_path": filter_vcf_path})


    def add_sg_variant_compare_effect(self, compare_id,file_path):
        """
        变异位点差异功能统计表
        """
        compare_id = self.check_objectid(compare_id)
        self.check_exists(file_path)
        data_list = []
        effect_type = [] # 后续在变异位点差异详情表中需要列举该类型，所以需要更新到主表。
        with open(file_path, 'r') as r:
            lines = r.readlines()
            for line in lines:
                item = line.strip().split("\t")
                insert_data = {
                "effect_type": item[0],
                "snp_num": int(item[1]),
                "indel_num": int(item[2]),
                "total_num": int(item[3]),
                "compare_id": compare_id
                }
                effect_type.append(item[0])
                data_list.append(insert_data)
        if len(data_list) == 0:
            self.bind_object.logger.info("{}文件为空！".format(file_path))
        else:
            self.col_insert_data("sg_variant_compare_effect", data_list)
            self.update_db_record("sg_variant_compare", {"_id": compare_id}, {"effect_type": effect_type})

    def add_variant_compare_effect_bar(self, compare_id, task_id, file_path, name = "all"):
        """
        sg_bar
        name 参数用于选择生成的snp，indel还是all。
        """
        compare_id = self.check_objectid(compare_id)
        self.check_exists(file_path)
        x_category = []
        value = []
        with open(file_path, 'r') as r:
            lines = r.readlines()
            for line in lines:
                item = line.strip().split("\t")
                x_category.append(item[0])
                if name == "snp":
                    value.append(int(item[1]))
                elif name == "indel":
                    value.append(int(item[2]))
                elif name == "all":
                    value.append(int(item[3]))
                else:
                    self.bind_object.logger.info("类型输入不正确")
        bar_id = self.add_sg_bar(compare_id, task_id, x_category, "variant_compare_effect", name)
        self.add_sg_bar_detail(bar_id, name, value)

    def add_sg_bar(self, feature_id, task_id, categories, location, name = "all"):
        """
        sg_bar
        """
        insert_data = {
            "origin_id": feature_id,
            "task_id": task_id,
            "location": location,
            "type": 1,
            "categories": categories,
            "name": name,
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

    def sg_variant_compare_impact(self, compare_id, file_path,):
        """
        变异位点差异功效统计表
        """
        compare_id = self.check_objectid(compare_id)
        self.check_exists(file_path)
        data_list = []
        with open(file_path, 'r') as r:
            lines = r.readlines()
            for line in lines:
                item = line.strip().split("\t")
                insert_data = {
                    "compare_id": compare_id,
                    "impact_type": item[0],
                    "snp_num": int(item[1]),
                    "indel_num": int(item[2]),
                    "total_num": int(item[3])}
                data_list.append(insert_data)
            if len(data_list) == 0:
                self.bind_object.logger.info("{}文件为空！".format(file_path))
            else:
                self.col_insert_data("sg_variant_compare_impact", data_list)

    def sg_variant_compare_impact_bar(self, compare_id, task_id, file_path, name = "all"):
        compare_id = self.check_objectid(compare_id)
        self.check_exists(file_path)
        data_list = []
        x_category = []
        value = []
        with open(file_path, 'r') as r:
            lines = r.readlines()
            for line in lines:
                item = line.strip().split("\t")
                x_category.append(item[0])
                if name == "snp":
                    value.append(int(item[1]))
                elif name == "indel":
                    value.append(int(item[2]))
                elif name == "all":
                    value.append(int(item[3]))
                else:
                    self.bind_object.logger.info("类型输入不正确")
        bar_id = self.add_sg_bar(compare_id, task_id, x_category, "variant_compare_impact", name )
        self.add_sg_bar_detail(bar_id, name, value)

    def sg_varian_compare_detail(self, compare_id, file_path, name = "all"):
        """
        变异位点差异详情表。
        """
        compare_id = self.check_objectid(compare_id)
        header = ["chr", "pos", "type", "ref", "alt", "ann"]  # 存放表头的list
        a = " "  # 用于测试
        with open(file_path, 'r') as r:
            lines = r.readlines()
            data_list = []
            for line in lines:
                # data_list = []
                table = {}  # 存放表头与内容的对应关系。
                insert_data = {}
                item = line.strip().split("\t")
                if re.match("#.", item[0]):
                    for j in item[6:]:
                        header.append(j)  # 标题中添加可变的内容
                else:
                    for row in range(6, len(item)):
                        if name == "snp":
                            if item[2] == "INDEL":
                                pass
                            else:
                                if header[row] not in table.keys():
                                    table[header[row]] = item[row]
                            insert_data = {
                                "compare_id": compare_id,
                                "chr": item[0],
                                "pos": item[1],
                                "type": item[2],
                                "ref": item[3],
                                "alt": item[4],
                                "ann": item[5],
                                "variant_type": name,
                            }
                        elif name == "indel":
                            if item[2] == "SNP":
                                pass
                            else:
                                if header[row] not in table.keys():
                                    table[header[row]] = item[row]
                            insert_data = {
                                "compare_id": compare_id,
                                "chr": item[0],
                                "pos": item[1],
                                "type": item[2],
                                "ref": item[3],
                                "alt": item[4],
                                "ann": item[5],
                                "variant_type": name,
                            }
                        else:
                            if header[row] not in table.keys():
                                table[header[row]] = item[row]
                            insert_data = {
                                "compare_id": compare_id,
                                "chr": item[0],
                                "pos": item[1],
                                "type": item[2],
                                "ref": item[3],
                                "alt": item[4],
                                "ann": item[5],
                                "variant_type": name,
                            }
                    for key1 in table.keys():
                            insert_data[key1] = table[key1]
                    data_list.append(insert_data)
            if len(data_list) == 0:
                self.bind_object.logger.info("{}文件为空！".format(file_path))
            else:
                print len(data_list)
                self.col_insert_data("sg_variant_compare_detail", data_list)
        self.update_db_record("sg_variant_compare", {"_id": compare_id}, {"header": header})

if __name__ == "__main__":
    a = VariantCompare(None)
    task_id = "tsg_32120"
    project_sn = "evolution_test"
    compare_id = "5b03bf9fa4e1af1482e33207"
    # a.sg_variant_compare("evolution_test", "")
    # a.sg_varian_compare_detail("5b03bf9fa4e1af1482e33520", "/mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/GeneticEvolution/new1/pop.table")
    # a.add_variant_compare_effect_bar("5b03bf9fa4e1af1482e33520", "evolution_test", "/mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/GeneticEvolution/new1/eff.type","all")
    # a.sg_variant_compare_impact_bar("5b03bf9fa4e1af1482e33520", "evolution_test", "/mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/GeneticEvolution/new1/function.type","snp")
    # a.sg_varian_compare_detail("5be942c0a4e1af0d6db35693", "/mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/GeneticEvolution/new1/pop.table", "all")
    a.sg_varian_compare_detail("5c0a6c5da4e1af033ea453cc","/mnt/ilustre/users/sanger-dev/workspace/20181210/VariantCompare_tsg_32120_1210103229361505_4103/VariantCompare/output/pop.table", "all")
    # a.add_sg_snp_indel_compare_stat(task_id, b1, "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/WGS/gatk_test/hongdong/sample/Ann.stat")
    # a.add_sg_snp_indel_compare_eff_stat(task_id, b1, "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/WGS/gatk_test/hongdong/sample/Eff.stat")
    # a.add_sg_snp_indel_compare_eff_stat(task_id, b2, "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/WGS/gatk_test/hongdong/sample/Eff.stat", "indel")
    # a.add_sg_snp_indel_compare_detail(b1, "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/WGS/gatk_test/hongdong/sample/pop.variant")
    # a.add_sg_snp_indel_compare_detail(compare_id, "/mnt/ilustre/users/sanger-dev/workspace/20180523/Single_tsg_29900_0523141120277624_4024/DoubleGroupCompare/output/diff.variant", "snp", "two_group")
    # compare_id = "5b0648e9a4e1af1f77f28b06"
    # results_path = "/mnt/ilustre/users/sanger-dev/workspace/20180524/Single_tsg_29900_0524130857538605_2512/DoubleGroupCompare/output/diff.variant.chr_1"
    # a.add_sg_manhattan(compare_id, results_path, types="indel")
    # results_path = "/mnt/ilustre/users/sanger-dev/workspace/ 20180414/Single_double_group_compare/DoubleGroupCompare/diff/win.stat"
    # compare_id = "5b03bf9fa4e1af1482e33207"
    # results_path = "/mnt/ilustre/users/sanger-dev/workspace/20180523/Single_tsg_29900_0523175002121620_8874/SingleGroupCompare/output/win.stat"
    # a.add_sg_distribution(compare_id, results_path, types="snp")
