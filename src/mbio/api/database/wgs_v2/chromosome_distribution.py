# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# modified 20180409

from api_base import ApiBase
import os
import re
import datetime
from collections import defaultdict


class ChromosomeDistribution(ApiBase):
    def __init__(self, bind_object):
        """
        染色体分布导表
        """
        super(ChromosomeDistribution, self).__init__(bind_object)
        self._project_type = "dna_wgs_v2"

    def add_sg_marker_distribution(self, project_sn, task_id, params=None):
        main_id = self.add_main_table("sg_marker_distribution", task_id, project_sn, params,
                                      "origin_chromosome_distribution", "染色体分布主表")
        return main_id

    def add_sg_area_detail(self, origin_id, file_path, project_sn, task_id, xmax):
        """
        添加基因组覆盖度分布细节表
        :param origin_id: 主表sg_area表id
        :param file_path: win_500000000_result.txt
        :param project_sn:
        :param task_id
        :param
        """
        origin_id = self.check_objectid(origin_id)
        area_id = self.sg_area(project_sn, task_id, origin_id, "")
        self.check_exists(file_path)
        data = defaultdict(list)
        legend = self.get_area_legend(file_path)
        if not legend:
            self.bind_object.logger.info("文件{}为空，不进行导表！".format(file_path))
            return
        with open(file_path, "r") as f:
            for line in f:
                if not re.match('#.*', line):
                    item = line.strip().split("\t")
                    if item[0] in legend:
                        data[item[0]] = [int(i) for i in item[1].split(',')]   # 数组字符元素转为int型
        self.sg_area_detail(area_id, legend, data)
        self.update_db_record("sg_area", {"origin_id": origin_id}, {"location": 'chromosome_distribution_area',
                                                                    "x_max": xmax})

    def get_area_legend(self, area_path):
        """
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
            legend = sca_list[:20]
        return legend

    def add_sg_distribution_detail(self, origin_id, file_path, project_sn, task_id, chr_file):
        """
        染色体分布图导表
        results_path: win.stat.xls
        """
        origin_id = self.check_objectid(origin_id)  # 检查id是否是OBjectID
        self.check_exists(file_path)  # 检查文件是否存在
        try:
            end = self.get_chr_max(chr_file)
        except:
            end = 100000
        main_data = [{
            "origin_id": origin_id,
            "name": "染色体分布图",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "location": 'chromosome_distribution_area',
            "type": 1,
            "project_sn": project_sn,
            "task_id": task_id
        }]
        main_id = self.col_insert_data("sg_distribution", main_data)
        legend = self.get_area_legend(file_path)
        insert_data = []
        win_data = defaultdict(list)
        with open(file_path, "r")as r:
            r.next()
            for line in r:
                temp = line.strip().split('\t')
                if temp[0] in legend:
                    win_data[temp[0]] = [int(i) for i in temp[1].split(',')]  # 数组字符元素转为int型

        for chr_ in legend:
            data = {
                "distribution_id": main_id,
                "name": chr_,
                "data_list": win_data[chr_]
            }
            insert_data.append(data)
            # if max(win_data[chr_]) > end:
            #     end = max(win_data[chr_])
        if len(insert_data) == 0:
            self.bind_object.logger.info("{}文件为空！".format(file_path))
        else:
            self.col_insert_data("sg_distribution_detail", insert_data)
            self.update_db_record("sg_distribution", {"_id": main_id},
                                  {"chr_list": legend, "start": 0, "end": end})

    def get_chr_max(self, file_path):
        chr_xmax = 0
        sca_xmax = 0
        with open(file_path, 'r') as r:
            for line in r:
                temp = line.strip().split('\t')
                if re.match('chr.*', temp[0]):
                    if int(temp[1]) > chr_xmax:
                        chr_xmax = int(temp[1])
                elif re.match('sca.*', temp[0]):
                    if int(temp[1]) > sca_xmax:
                        sca_xmax = int(temp[1])
        result = chr_xmax
        if chr_xmax == 0:
            result = sca_xmax
        return result



if __name__ == "__main__":
    a = ChromosomeDistribution(None)
