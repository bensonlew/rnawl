# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 20180916

from api_base import ApiBase
import datetime
from collections import defaultdict
from bson.objectid import ObjectId
import re
import gzip
from types import StringTypes


class LdAnalysis(ApiBase):
    def __init__(self, bind_object):
        """
        连锁不平衡导表
        """
        super(LdAnalysis, self).__init__(bind_object)
        self._project_type = "dna_evolution"

    def sg_ld(self, task_id, project_sn, params):
        """
        变异位点比较分析主表
        :return:
        """
        main_id = self.add_main_table("sg_ld", task_id, project_sn, params, "origin_ld_analysis", "连锁不平衡接口",
                                      "开始进行连锁不平衡插图")
        self.update_db_record("sg_ld", {"_id": main_id}, {"main_id": main_id, "status": "end"})
        return main_id

    def update_ld_analysis(self, graph_path, update_id):
        update_id = self.check_objectid(update_id)
        self.update_db_record("sg_ld", {"_id": update_id}, {"graph_path": graph_path
                                                                    })

    def sg_ld_detail(self, ld_id, file_path):
        ld_id = self.check_objectid(ld_id)
        with open(file_path)as r:
            header = []
            dict1 = {}
            data_list1 = []
            data_list2 = []
            data_list3 = []
            lines = r.readlines()
            for i in lines[1:]:
                list1 = []
                group_name = i.strip().split("\t")[0]
                header.append(group_name)
                r2_8_index = 0  # 第一列的数字刚开始赋值为0。
                r2_8_data = 3  # 这三个字都赋值为-3，保证肯定不会取到这三个值。
                r2_1_index = 0
                r2_1_data = 3
                r2_5_index = 0
                r2_5_data = 3
                max1 = 0
                half_life = 0
                with gzip.open(i.strip().split("\t")[1], "r") as f:
                    ld_lines = f.readlines()
                    for k in ld_lines[1:]:
                        if float(k.split("\t")[1]) > max1:
                            max1 = float(k.split("\t")[1])
                            half_life = round(max1 / float(2), 4)
                    for j in ld_lines[1:]:
                        if float(j.split("\t")[1])-float(0.8) >= 0:
                            if float(r2_8_data)-float(0.8) > float(j.split("\t")[1]) - float(0.8):
                                r2_8_data = j.split("\t")[1]
                                r2_8_index = j.split("\t")[0]
                        if float(j.split("\t")[1]) - float(0.1) >= 0:
                            if float(r2_1_data) - float(0.1) > float(j.split("\t")[1]) - float(0.1):
                                r2_1_data = j.split("\t")[1]
                                r2_1_index = j.split("\t")[0]
                        if float(j.split("\t")[1]) - half_life >= 0:
                            if float(r2_5_data) - half_life > float(j.split("\t")[1]) - half_life:
                                r2_5_data = j.split("\t")[1]
                                r2_5_index = j.split("\t")[0]
                    if r2_8_index == 0:
                        r2_8_index = "/ "
                    list1.append(r2_8_index)
                    if r2_1_index == 0:
                        r2_1_index = "/ "
                    list1.append(r2_1_index)
                    if r2_5_index == 0:
                        r2_5_index = "/ "
                    list1.append(r2_5_index)
                dict1[i.strip().split("\t")[0]] = list1
            table_r28 = {}
            table_r21 = {}
            table_r25 = {}
            for x in header:
                table_r28[x] = dict1[x][0]
                table_r21[x] = dict1[x][1]
                table_r25[x] = dict1[x][2]
            # for n in ["r2_8", "r2_1", "r2_5"]:
            for n in ["R2 > 0.8", "R2 > 0.1", "R2 half-life"]:
                insert_data = {"pops": n, "ld_id": ld_id}
                if n == "R2 > 0.8":
                    for key1 in table_r28.keys():
                        insert_data[key1] = table_r28[key1]
                    data_list1.append(insert_data)
                    if len(data_list1) == 0:
                        self.bind_object.logger.info("r2_8{}文件为空！".format(file_path))
                    else:
                        self.col_insert_data("sg_ld_detail", data_list1)
                elif n == "R2 > 0.1":
                    for key2 in table_r21.keys():
                        insert_data[key2] = table_r21[key2]
                    data_list2.append(insert_data)
                    if len(data_list2) == 0:
                        self.bind_object.logger.info("r2_1{}文件为空！".format(file_path))
                    else:
                        self.col_insert_data("sg_ld_detail", data_list2)
                elif n == "R2 half-life":
                    for key3 in table_r25.keys():
                        insert_data[key3] = table_r25[key3]
                    data_list3.append(insert_data)
                    if len(data_list3) == 0:
                        self.bind_object.logger.info("r2_1{}文件为空！".format(file_path))
                    else:
                        self.col_insert_data("sg_ld_detail", data_list3)
            self.update_db_record("sg_ld", {"_id": ld_id}, {"header": header})

    # def add_ld_curve(self, file_path, ld_curve_id, task_id):
    #     ld_curve_id = self.check_objectid(ld_curve_id)
    #     gro_dict = {}
    #     # gro_list = []
    #     with open(file_path)as r:
    #         lines = r.readlines()
    #         for i in lines[1:]:
    #             if i.strip().split("\t")[0] not in gro_dict.keys():
    #                 gro_dict[i.strip().split("\t")[0]] = []
    #                 # with gzip.open(i.strip().split("\t")[1], "r") as f:
    #                 with open(i.strip().split("\t")[1], "r") as f:
    #                     ld_lines = f.readlines()
    #                     for j in ld_lines[1:]:
    #                         # coordinate = [j.strip().split()[0], j.strip().split()[1]]
    #                         coordinate = float(j.strip().split()[2])
    #                         gro_dict[i.strip().split("\t")[0]].append(coordinate)
    #     curve_id = self.sg_curve(task_id, ld_curve_id, "LD-decay Distribution", "", 1, "ld_decay_curve")  # 这个数字需要根据前端修改
    #     for key1 in gro_dict.keys():
    #         self.sg_curve_detail(curve_id, key1, gro_dict[key1])

    # def sg_curve(self, task_id, origin_id, name, categories, types, location=None, other_attr=None, ext_name=None):
    #     """
    #     用到导曲线图的数据
    #     :param task_id:
    #     :param origin_id: 上一级主表id，如sg_specimen_qc_id,sg_mapping_detail_id等
    #     :param name: 图的中文名称
    #     :param categories: x轴
    #     :param types: 类型，普通曲线
    #     :param location: 位置，对应的图英文名称，用于定位当前关联ID的哪个位置，通常用于解决一个结果表画相同的类型的图
    #     :param other_attr: 用于记录一些标记
    #     :param ext_name: 用于标题中显示名称的，当type为sample，改名时标题也需要改名：
    #     example：{"type":"sample","title":"样本原名"，"ext_title":"_1"}
    #     :return:
    #     """
    #     if not isinstance(origin_id, ObjectId):
    #         if isinstance(origin_id, StringTypes):
    #             origin_id = ObjectId(origin_id)
    #         else:
    #             raise Exception("origin_id必须为ObjectId对象或其对应的字符串!")
    #     data_list = []
    #     insert_data = {
    #         "task_id": task_id,
    #         "origin_id": origin_id,
    #         "name": name,
    #         "categories": categories,
    #         "type": types,
    #         "location": location if location else "",
    #         "other_attr": other_attr if other_attr else "",
    #         "ext_name": ext_name if ext_name else "",
    #         "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    #     }
    #     data_list.append(insert_data)
    #     main_id = self.col_insert_data("sg_curve", data_list)
    #     return main_id
    #
    # def sg_curve_detail(self, curve_id, name, value, is_show_log="true"):
    #     """
    #     添加曲线图细节表
    #     :param curve_id: 曲线图的主表id
    #     :param name:  样本名称
    #     :param value:  值不能是字符串，值的顺序按照主表的categories里面对应的顺序
    #     :param is_show_log:  true  or false 用于决定是否显示插入成功的log日志，true为显示
    #     :return:
    #     """
    #     data_list = []
    #     insert_data = {
    #         "curve_id": curve_id,
    #         "name": name,
    #         "value": value
    #     }
    #     data_list.append(insert_data)
    #     self.col_insert_data("sg_curve_detail", data_list, is_show_log)

if __name__ == "__main__":
    a = LdAnalysis(None)
    task_id = "evolution_test"
    project_sn = "evolution_test"
    compare_id = "5b03bf9fa4e1af1482e33207"
    a.sg_ld(project_sn, "我是参数")
    # a.sg_ld_detail("5b03bf9fa4e1af1482e33207", "/mnt/ilustre/users/sanger-dev/workspace/20180917/Single_test_ld_delay20180917153221/LdAnalysis/gro_list")
    # a.add_lg_curve("/mnt/ilustre/users/sanger-dev/workspace/20180917/Single_test_ld_delay20180917153221/LdAnalysis/gro_list", "5b03bf9fa4e1af1482e33207")
#     a.sg_variant_compare("woshishi", "canshu shi shadongxi ")
#     # a.sg_varian_compare_detail("5b03bf9fa4e1af1482e33520", "/mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/GeneticEvolution/new1/pop.table")
#     a.add_variant_compare_effect_bar("5b03bf9fa4e1af1482e33520", "evolution_test", "/mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/GeneticEvolution/new1/eff.type","ALL")
#     a.sg_variant_compare_impact_bar("5b03bf9fa4e1af1482e33520", "evolution_test", "/mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/GeneticEvolution/new1/function.type","SNP")
#     a.sg_varian_compare_detail("5b03bf9fa4e1af1482e33520", "/mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/GeneticEvolution/new1/pop.table", "SNP")
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
