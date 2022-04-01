# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 20190304

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
        self._project_type = "dna_wgs_v2"

    def sg_variant_compare(self, task_id, project_sn, params):
        """
        变异位点比较分析主表
        :return:
        """
        main_id = self.add_main_table("sg_variant_compare", task_id, project_sn, params, "origin_variant_compare",
                                      "比较分析接口导表", "开始进行有参变异位点比较插图")
        self.update_db_record("sg_variant_compare", {"_id": main_id}, {"main_id": main_id, "status": "end"})
        return main_id

    def update_variant_compare(self, filter_vcf_path, update_id):
        update_id = self.check_objectid(update_id)
        self.update_db_record("sg_variant_compare", {"_id": update_id}, {"vcf_path": filter_vcf_path})

    def add_sg_variant_compare_stat(self, compare_id, file_path, name):
        """
        这里传入的路径为过滤后的vcf格式
        :param compare_id:
        :param file_path:
        :return:
        """
        compare_id = self.check_objectid(compare_id)
        chr_snp = {}  # 用于存放snp个数
        chr_indel = {}  # 用于存放indel格式
        index_snp = 0  # 用于表示snp的个数的指针
        index_indel = 0  # 用于表示indel的个数的指针
        with open(file_path) as f:
            lines = f.readlines()
            for i in range(len(lines) - 1):  # 这里最后一行做单独判断
                if re.match("#", lines[i]):
                    pass
                else:
                    variant = lines[i].strip().split("\t")
                    gene_list = (str(variant[3]) + "," + str(variant[4])).strip().split(",")
                    variant_type = "SNP"
                    for gene in gene_list:
                        if len(gene) > 1:
                            variant_type = "INDEL"
                            break
                        else:
                            pass
                    if variant_type == "INDEL":
                        if lines[i].strip().split("\t")[0] == lines[i + 1].strip().split("\t")[0]:
                            index_indel += 1
                        else:
                            if not lines[i].strip().split("\t")[0] in chr_indel.keys():
                                chr_indel[lines[i].strip().split("\t")[0]] = index_indel
                                index_indel = 0
                    else:
                        if lines[i].strip().split("\t")[0] == lines[i + 1].strip().split("\t")[0]:
                            index_snp += 1
                        else:
                            if not lines[i].strip().split("\t")[0] in chr_snp.keys():
                                chr_snp[lines[i].strip().split("\t")[0]] = index_snp
                                index_snp = 0
            variant = lines[len(lines) - 1].strip().split("\t")
            gene_list = (str(variant[3]) + "," + str(variant[4])).strip().split(",")
            variant_type = "SNP"
            for gene in gene_list:
                if len(gene) > 1:
                    variant_type = "INDEL"
                    break
                else:
                    pass
            if variant_type == "INDEL":
                if not lines[len(lines) - 1].strip().split("\t")[0] in chr_indel.keys():
                    chr_indel[lines[len(lines) - 1].strip().split("\t")[0]] = index_indel + 1  # 如果和倒数第二个不一样，前面
                    chr_snp[lines[len(lines) - 1].strip().split("\t")[0]] = index_snp
                    # 已经结算，
                    # 这就就是1，如果一样，前面还没有结算，可以继续加1
            else:
                if not lines[len(lines) - 1].strip().split("\t")[0] in chr_snp.keys():
                    chr_snp[lines[len(lines) - 1].strip().split("\t")[0]] = index_snp + 1
                    chr_indel[lines[len(lines) - 1].strip().split("\t")[0]] = index_indel
        chr_all = []
        data_list = []
        data_list1 = []
        for keys in chr_snp.keys():
            if keys not in chr_all:
                chr_all.append(keys)
        for keys1 in chr_indel.keys():
            if keys1 not in chr_all:
                chr_all.append(keys1)
        for items in chr_all:
            if items not in chr_indel.keys():
                chr_indel[items] = 0
            elif items not in chr_snp.keys():
                chr_snp[items] = 0
            insert_data = {
                "chr_id": items,
                "snp_num": chr_snp[items],
                "indel_num": chr_indel[items],
                "total_num": chr_snp[items] + chr_indel[items],
                "name": name,
                "compare_id": compare_id,
            }
            """
            用于结果的排序"""
            if items.startswith("chr"):
                data_list.append(insert_data)
            elif items.startswith("sca"):
                data_list1.append(insert_data)
            else:
                raise Exception("请检查染色体格式")

        if len(data_list) == 0 and len(data_list1) == 0:
            self.bind_object.logger.info("{}文件为空！".format(file_path))
        else:
            # data_list.sort(key=lambda k: int(re.findall("\d+", k["chr_id"])[0]))
            # data_list1.sort(key=lambda k: int(re.findall("\d+", k["chr_id"])[0]))
            data_list.sort(key=lambda k: k["chr_id"])
            data_list1.sort(key = lambda k: k["chr_id"])
            new_list = data_list + data_list1
            self.col_insert_data("sg_variant_compare_stat", new_list)

    def add_sg_variant_compare_stat_v2(self, compare_id, file_path, name):
        compare_id = self.check_objectid(compare_id)
        data_list = []
        with open(file_path, 'r') as r:
            for line in r:
                if re.match('#', line):
                    pass
                else:
                    temp = line.strip().split('\t')
                    insert_data = {
                        "chr_id": temp[0],
                        "snp_num": int(temp[1]),
                        "indel_num": int(temp[2]),
                        "total_num": int(temp[3]),
                        "name": name,
                        "compare_id": compare_id,
                    }
                    data_list.append(insert_data)
        if len(data_list) == 0:
            self.bind_object.logger.info("{}文件为空！".format(file_path))
        else:
            self.col_insert_data("sg_variant_compare_stat", data_list)

    def add_sg_variant_compare_effect(self, compare_id,file_path, name):
        """
        变异位点差异功能统计表
        """
        compare_id = self.check_objectid(compare_id)
        self.check_exists(file_path)
        data_list = []
        effect_type = [] # 后续在变异位点差异详情表中需要列举该类型，所以需要更新到主表。
        with open(file_path, 'r') as r:
            lines = r.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                "effect_type": item[0],
                "snp_num": int(item[1]),
                "indel_num": int(item[2]),
                "total_num": int(item[3]),
                "compare_id": compare_id,
                    "name": name,
                }
                effect_type.append(item[0])
                data_list.append(insert_data)
        if len(data_list) == 0:
            self.bind_object.logger.info("{}文件为空！".format(file_path))
        else:
            self.col_insert_data("sg_variant_compare_effect", data_list)
            self.update_db_record("sg_variant_compare", {"_id": compare_id}, {"effect_type": effect_type})

    def add_variant_compare_effect_bar(self, compare_id, task_id, file_path, subname, name = "all"):  # 这里的subname用于二级弹窗
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
            for line in lines[1:]:
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
        bar_id = self.add_sg_bar(compare_id, task_id, x_category, "variant_compare_impact", name)
        self.update_db_record("sg_bar", {"_id": bar_id}, {"name": subname})
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
            "variant_type": name,
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

    def add_sg_variant_compare_func(self, compare_id, file_path, name):
        """
        变异位点差异功效统计表
        """
        compare_id = self.check_objectid(compare_id)
        self.check_exists(file_path)
        data_list = []
        func_type = []
        with open(file_path, 'r') as r:
            lines = r.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "compare_id": compare_id,
                    "func_type": item[0],
                    "snp_num": int(item[1]),
                    "indel_num": int(item[2]),
                    "total_num": int(item[3]),
                    "name": name
                }
                func_type.append(item[0])
                data_list.append(insert_data)
            if len(data_list) == 0:
                self.bind_object.logger.info("{}文件为空！".format(file_path))
            else:
                self.col_insert_data("sg_variant_compare_func", data_list)
                self.update_db_record("sg_variant_compare", {"_id": compare_id}, {"func_type": func_type})

    def sg_variant_compare_impact_bar(self, compare_id, task_id, file_path, subname, name = "all"):
        compare_id = self.check_objectid(compare_id)
        self.check_exists(file_path)
        data_list = []
        x_category = []
        value = []
        with open(file_path, 'r') as r:
            lines = r.readlines()
            for line in lines[1:]:
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
        bar_id = self.add_sg_bar(compare_id, task_id, x_category, "variant_compare_func", name )
        self.update_db_record("sg_bar", {"_id": bar_id}, {"name": subname})
        self.add_sg_bar_detail(bar_id, name, value)

    def sg_varian_compare_detail(self, compare_id, file_path, name):
        """
        变异位点差异详情表。
        """
        compare_id = self.check_objectid(compare_id)
        header = []
        header_dict = {}
        with open(file_path, 'r') as r:
            lines = r.readlines()
            data_list = []
            for line in lines:
                table = {}  # 存放表头与内容的对应关系。
                insert_data = {}
                item = line.strip().split("\t")
                if re.match("#.", item[0]):
                    for j in item[8:]:
                        header.append(j)  # 标题中添加可变的内容
                else:
                    for row in range(8, len(item)):
                        if header[row-8] not in table.keys():
                            table[header[row-8]] = item[row]
                    insert_data = {
                        "compare_id": compare_id,
                        "chr": item[0],
                        "pos": int(item[1]),
                        "snp/indel_id": item[2] if not item[2] != "" else str(item[0])+"_"+ str(item[1]),
                        "ref": item[3],
                        "alt": item[4],
                        "type": item[5],
                        "name": name,
                        "annotation": item[6],
                        "region": item[7],
                            }
                    for key1 in table.keys():
                        if re.match(".*_genotype", key1):
                            header_dict[key1] = key1.strip().split("_")[0] + " " + "Genotype"
                        elif re.match(".*_allele_depth", key1):
                            header_dict[key1] = key1.strip().split("_")[0] + " " + "Allele Depth"
                        elif re.match(".*_alle_fre", key1):
                            header_dict[key1] = key1.strip().split("_")[0] + " " + "Allele Frequency"
                        else:
                            print "________没有这种情况_____________________"
                        insert_data[key1] = table[key1]
                    data_list.append(insert_data)
            if len(data_list) == 0:
                self.bind_object.logger.info("{}文件为空！".format(file_path))
            else:
                self.col_insert_data("sg_variant_compare_detail", data_list)
        if name == "variant_compare":  # 只有当时mutiply的时候，name才有可能是variant_compare
            self.update_db_record("sg_variant_compare", {"_id": compare_id}, {"header": header})
            self.update_db_record("sg_variant_compare", {"_id": compare_id}, {"header_dict": header_dict})

            """
            添加一个header_dict的字段并更新到主表"""

if __name__ == "__main__":
    a = VariantCompare(None)
    task_id = "wgs_v2_new"
    project_sn = "wgs_v2"
    # a.sg_variant_compare(task_id,project_sn, "wgs_v2 test")
    name = "variant_compare"
    subname = "16S15_vs_16S456_1_compare"
    compare_id = "5c860c5117b2bf7c7a2db603"
    a.add_sg_variant_compare_stat(compare_id, "/mnt/ilustre/users/sanger-dev/i-sanger_workspace/20190718/"
                                              "VariantCompare_i-sanger_190269_0718084910770247_5927/output/variant_compare/W305_vs_S1B_compare.filter.vcf", name)
    # a.add_sg_variant_compare_stat(compare_id, "/mnt/ilustre/users/sanger-dev/workspace/20190418/VariantCompare_sanger_85433_0418174033526279_4095/output/variant_compare/GC_bulk_vs_ZH30_compare.filter.vcf", name)
    a.add_sg_variant_compare_effect(compare_id,"/mnt/ilustre/users/sanger-dev/i-sanger_workspace/20190718/"
                                              "VariantCompare_i-sanger_190269_0718084910770247_5927/output/variant_compare/W305_vs_S1B_compare.eff", name)
    a.add_sg_variant_compare_func(compare_id, "/mnt/ilustre/users/sanger-dev/i-sanger_workspace/20190718/"
                                              "VariantCompare_i-sanger_190269_0718084910770247_5927/output/variant_compare/W305_vs_S1B_compare.func", name)
    a.sg_varian_compare_detail(compare_id,"/mnt/ilustre/users/sanger-dev/i-sanger_workspace/20190718/"
                                              "VariantCompare_i-sanger_190269_0718084910770247_5927/output/variant_compare/W305_vs_S1B_compare.detail", name)
    a.sg_variant_compare_impact_bar(compare_id, task_id,"/mnt/ilustre/users/sanger-dev/i-sanger_workspace/20190718/"
                                              "VariantCompare_i-sanger_190269_0718084910770247_5927/output/variant_compare/W305_vs_S1B_compare.eff", subname, "all")
    a.sg_variant_compare_impact_bar(compare_id, task_id, "/mnt/ilustre/users/sanger-dev/i-sanger_workspace/20190718/"
                                              "VariantCompare_i-sanger_190269_0718084910770247_5927/output/variant_compare/W305_vs_S1B_compare.func", subname, "all")

