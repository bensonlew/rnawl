# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 20180827

from api_base import ApiBase
import datetime
from collections import defaultdict
import re
import os
import math
import json


class ChromosomeWindow(ApiBase):
    def __init__(self, bind_object):
        """
        基因型分布导表

        """
        super(ChromosomeWindow, self).__init__(bind_object)
        self._project_type = "dna_evolution"
        # self.project_sn = self.bind_object.sheet.project_sn
        # self.task_id = self.bind_object.sheet.id
        # self.project_sn = "188_5bbc1f5fd79d2"
        # self.task_id = "tsg_32120"

    def add_sg_chromosome_window(self, project_sn, params, task_id, compare_id):
        """
        基因组比对基因组覆盖分布图导表
        """
        compare_id = self.check_objectid(compare_id)
        # main_id = self.add_main_table("sg_chromosome_window", task_id, project_sn, params, "",
        #                                     "chromosome_window导表",
        #                                     "开始进行chromosome_window插图导表")
        data_list = []
        insert_data = {
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "desc": " 基因组比对基因组覆盖分布图导表",
            "task_id": task_id,
            "project_sn": project_sn,
            "name": "origin_chromosome_window",
            "params": json.dumps(params, sort_keys=True, separators=(',', ':'))
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_chromosome_window", data_list)
        self.update_db_record("sg_chromosome_window", {"_id": main_id}, {"main_id": main_id, "compare_id": compare_id})
        return main_id

    def add_sg_distribution(self, path_dir, main_id, location):
        main_id = self.check_objectid(main_id)
        """
        location 用于标记是snp，还是indel 还是 all。
        :param path_file:
        :param main_id:
        :param location:
        :return:
        """
        file_path = ""
        if location == "snp":
            file_path = os.path.join(path_dir, "data_snp")
        elif location == "indel":
            file_path = os.path.join(path_dir, "data_indel")
        elif location == "all":
            file_path = os.path.join(path_dir, "data_all")
        insert_data = []
        main_data = [{
            "origin_id": main_id,
            "name": "染色体分布图",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "location": "variant_distribution",
            "variant_type": location,
            "type": 1
        }]
        distribution_id = self.col_insert_data("sg_distribution", main_data)
        with open(file_path) as r:
            lines = r.readlines()
            final_list = lines[0].strip().split(",")
            win_data = json.loads(lines[1])
            start = lines[2]
            end = lines[3]
            for chr in final_list:
                data = {
                    "distribution_id": distribution_id,
                    "name": chr,
                    "data_list": win_data[chr]
                }
                insert_data.append(data)
            if len(insert_data) == 0:
                self.bind_object.logger.info("{}文件为空！".format(path_dir))
            else:
                self.col_insert_data("sg_distribution_detail", insert_data)
                self.update_db_record("sg_distribution", {"_id": distribution_id},
                                      {"chr_list": final_list, "start": start, "end": end})


if __name__ == "__main__":
    a = ChromosomeWindow(None)
    task_id = "tsg_32120"
    project_sn = "188_5bbc1f5fd79d2"
    compare_id = "5bbef71da4e1af72e2c284ea"
    path_dir = '/mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/GeneticEvolution/try'
    # a.add_sg_chromosome_window(project_sn,"", task_id, compare_id)
    a.add_sg_distribution(path_dir, "5bbf0c4ba4e1af2cb87ea8ba", "all")
    a.add_sg_distribution(path_dir, "5bbf0c4ba4e1af2cb87ea8ba", "snp")
    a.add_sg_distribution(path_dir, "5bbf0c4ba4e1af2cb87ea8ba", "indel")
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
