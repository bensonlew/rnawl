# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 20180827

from api_base import ApiBase
import datetime
from collections import defaultdict
import re
import os
import math


class CoverageWindow(ApiBase):
    def __init__(self, bind_object):
        """
        基因组分布导表
        """
        super(CoverageWindow, self).__init__(bind_object)
        self._project_type = "dna_wgs_v2"
        # self.project_sn = self.bind_object.sheet.project_sn
        # self.task_id = self.bind_object.sheet.id
        # self.project_sn = "woshi project_sn"
        # self.task_id = "woshi task_id"
        self.name = ''

    def add_sg_coverage_window(self, params, mapping_id):
        """
        基因组比对基因组覆盖分布图导表
        """
        mapping_id = self.check_objectid(mapping_id)
        main_id = self.add_main_table("sg_coverage_window", self.bind_object.sheet.id,
                                      self.bind_object.sheet.project_sn, params, "origin_coverage_window",  "coverage_window导表",
                                      "开始进行coverage_wind插图导表")
        self.update_db_record("sg_coverage_window", {"_id": main_id}, {"main_id": main_id, "mapping_id": mapping_id,
                                                                       "status": "end"})
        return main_id

    def get_sample(self, path_dir, coverage_id):
        sample_list = []  # 用于存放样本的list，存放于主表。
        coverage_id = self.check_objectid(coverage_id)
        for f in os.listdir(path_dir):
            self.name = f.split(".")[0]
            sample_list.append(self.name)
            area_id = self.add_sg_area(self.name, coverage_id)
            area_path = os.path.join(path_dir, f)
            self.add_sg_area_detail(area_id, area_path, coverage_id)
        self.update_db_record("sg_coverage_window", {"_id": coverage_id}, {"sample_list": sample_list})

    def add_sg_area(self, name, origin_id):
        """
        添加基因组覆盖度分布主表
        :param origin_id: 主表sg_coverage_window表id
        :param name: 样本名称
        """
        main1_id = self.sg_area(self.bind_object.sheet.project_sn, self.bind_object.sheet.id, origin_id, name)
        return main1_id

    def add_sg_area_detail(self, origin_id, area_path, coverage_id):
        """
        添加基因组覆盖度分布细节表
        :param coverage_id: 主表sg_area表id
        :param area_path: Lands.coverage
        :param
        """
        coverage_id = self.check_objectid(coverage_id)
        self.check_exists(area_path)
        data = {}
        legend = self.get_area_legend(area_path)
        with open(area_path, "r") as f:
            for line in f:
                item = line.strip().split("\t")
                if item[0] in legend:
                    if item[0] not in data.keys():
                        data[item[0]] = []
                    value = round(math.log(float(item[2])) / math.log(2), 4)
                    data[item[0]].append(value)
        self.sg_area_detail(origin_id, legend, data)
        self.update_db_record("sg_coverage_window", {"_id": coverage_id}, {"chr_list": legend}, is_show_log='false')

    def get_area_legend(self, area_path):
        """s
        得到基因组覆盖度图coverage文件的legend
        若有chr，则为所有的chr；若没有chr，则为所有的sca
        """
        chr_list, sca_list = [], []
        with open(area_path, "r") as f:
            for line in f:
                item = line.strip().split("\t")
                if item[0].startswith("chr"):
                    if item[0] not in chr_list:
                        chr_list.append(item[0])
                else:
                    if item[0] not in sca_list:
                        sca_list.append(item[0])
        if chr_list:
            legend = chr_list
        else:
            legend = sca_list
        return legend


if __name__ == "__main__":
    a = CoverageWindow(None)
    task_id = " tsg_32120"
    project_sn = "evolution_test"
    compare_id = "5bac7343a4e1af5f8dc8ab99"
    path_dir = '/mnt/ilustre/users/sanger-dev/workspace/20180905/Single_tsanger_30729_CoverageWindows_0905174710486344/CoverageWindows/output/coverage_window_dir'
    a.add_sg_coverage_window("I am params", "5badbb2fa4e1af0d69bb97d3")
    a.get_sample(path_dir, "5b03bf9fa4e1af1482e33207", "5badbb2fa4e1af0d69bb97d3")
    # a.sg_varian_compare_detail("5b03bf9fa4e1af1482e33520", "/mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/GeneticEvolution/new1/pop.table")
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
