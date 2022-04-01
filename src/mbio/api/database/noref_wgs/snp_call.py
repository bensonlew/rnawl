# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 20181220

import os
import types
import datetime
from bson.son import SON
from bson.objectid import ObjectId
from biocluster.api.database.base import Base, report_check
from api_base import ApiBase


class SnpCallApi(ApiBase):
    """
    snp call 和 snp compare 都在这里
    """
    def __init__(self, bind_object):
        super(SnpCallApi, self).__init__(bind_object)
        self._project_type = "dna_noref_wgs"

    def add_sg_snp_call(self, task_id, params):
        """
        sg_lg
        add by hongdong@20180705
        :return:
        """
        main_id = self.add_main_table("sg_snp_call", task_id, "", params, "snp_call",
                                      "snp变异检测主表")
        self.update_db_record("sg_snp_call", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_snp_call_stat(self, path, call_id):
        """
        snp数据统计表导表
        :param path:
        :return:
        """
        call_id = self.check_objectid(call_id)
        with open(path) as f:
            lines = f.readlines()
            data_list = []
            for line in lines[1:]:
                iterm = line.strip().split()
                insert_data ={
                    "call_id": call_id,
                    "specimen_id": iterm[0],
                    "num":iterm[1],
                    "transition":iterm[2],
                    "transversion":iterm[3],
                    "ti_tv":iterm[4],
                    "hete_num":iterm[5],
                    "homo_num":iterm[6]
                }
                data_list.append(insert_data)
            if len(data_list) != 0:
                self.col_insert_data("sg_snp_call_stat", data_list)
            else:
                self.bind_object.logger.info("文件{}为空，不进行导表！".format(path))

    def add_snp_depth_curve(self, task_id, call_id, file_path, location, name):
        """
        snp质量评估的质量分布图/深度分布图
        """
        origin_id = self.check_objectid(call_id)
        self.check_exists(file_path)
        samples = []
        data = {}
        min, max = 0, 0
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t") # 分隔符注意
                print item
                sample_id = item[0]
                if sample_id not in samples:
                    samples.append(sample_id)
                    data[sample_id] = []
                if int(item[1]) < min:
                    min = int(item[1])
                if int(item[1]) > max:
                    max = int(item[1])
                data[sample_id].append(int(item[2]))
        categories = []
        for i in range(min, max+1):
            categories.append(str(i))
        curve_id = self.sg_curve(task_id, origin_id, name, categories, 1, location, '', '')
        for s in samples:
            sum = float(data[s][-1])
            percent_data = []
            for n in data[s]:
                percent = round(float(n) / sum, 4) * 100                 # 确认这里是否需要计算，以及计算的方式是什么
                percent_data.append(percent)
            self.sg_curve_detail(curve_id, s, percent_data)
            print "导入样本{}snp深度分布曲线图成功".format(s)

    def add_snp_density_bar(self, task_id, call_id, file_path):
        """
        """
        name = "Tag SNP Density Distribution"
        categories = []
        value = []
        with open(file_path, "r")as fr:
            lines = fr.readlines()
            for line in lines[1:]:
                tmp = line.strip().split("\t")
                categories.append(tmp[0])
                value.append(tmp[1])
        bar_id = self.sg_bar(task_id, call_id, name, categories, 1, "snp_density")
        self.sg_bar_detail(bar_id, "snp_density", value)


if __name__ == "__main__":
    data = SnpCallApi(None)
    task_id = "test"
    params = "无参导表"
    call_id = "5c1b6eb3a4e1af55edd74c95"
    # data.add_sg_snp_call(task_id, params)
    # data.add_sg_snp_call_stat("/mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/wgs_noref/api/snp.stat", "5c1b6eb3a4e1af55edd74c95")
    # data.add_snp_depth_curve(task_id, call_id, "/mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/wgs_noref/api/"
    #                                         "curve/depth_distribution", "snp深度分布图", "snp_depth_distribution")
    data.add_snp_density_bar(task_id, call_id, "/mnt/ilustre/users/sanger-dev/sg-users/liuwentian/wucan/06.genotype/consensus/tag_snp.txt")




