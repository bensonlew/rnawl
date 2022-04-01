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



class SnpCompare(ApiBase):
    """
    比较分析脚本

    """
    def __init__(self, bind_object):
        super(SnpCompare, self).__init__(bind_object)
        self._project_type = "dna_noref_wgs"

    def add_sg_snp_compare(self, task_id, params):
        """
        :return:
        """
        main_id = self.add_main_table("sg_snp_compare", task_id, "", params, "snp_compare",
                                      "snp变异比较分析主表")
        self.update_db_record("sg_snp_compare", {"_id": main_id}, {"main_id": main_id})
        return main_id
    def add_sg_snp_compare_stat(self, path, compare_id, name):
        """
        snp数据统计表导表
        :param path:
        :return:
        """
        compare_id = self.check_objectid(compare_id)
        with open(path) as f:
            lines = f.readlines()
            data_list = []
            for line in lines[1:]:
                iterm = line.strip().split("\t")
                insert_data =\
                {
                    "snp_num":iterm[0],
                    "average_depth":iterm[1],
                    "miss_ratio":iterm[2],
                    "compare_id":compare_id,
                    "name":name
                }
                data_list.append(insert_data)
            if len(data_list) != 0:
                self.col_insert_data("sg_snp_compare_stat", data_list)
            else:
                self.bind_object.logger.info("文件{}为空，不进行导表！".format(path))

    def add_snp_compare_detail(self, path, compare_id, name):
        header = ["SNP ID", "Consensus ID", "Pos", "Ref"]  # 存放表头的list
        compare_id= self.check_objectid(compare_id)
        with open(path, 'r') as r:
            lines = r.readlines()
            data_list = []
            for line in lines:
                table = {}  # 存放表头与内容的对应关系。
                item = line.strip().split("\t")
                if lines.index(line) ==0:
                    for j in item[5:-1]:
                        header.append(j)  # 标题中添加可变的内容
                else:
                    insert_data = {
                        "name": name,
                        "compare_id": compare_id,
                        "Consensus ID": item[0],
                        "Pos": item[1],
                        "SNP ID": item[2],
                        "Ref": item[3],
                    }
                    for row in range(4, len(item) -2): # 最后一列不要
                        if header[row] not in table.keys():
                            table[header[row]] = item[row+1]
                    for key1 in table.keys():
                        insert_data[key1] = table[key1]
                    data_list.append(insert_data)

            if len(data_list) == 0:
                self.bind_object.logger.info("{}文件为空！".format(path))
            else:
                self.col_insert_data("sg_snp_compare_detail", data_list)
        self.update_db_record("sg_snp_compare", {"_id": compare_id}, {"header": header})

    def update_snp_compare(self, filter_vcf_path, update_id):
        update_id = self.check_objectid(update_id)
        self.update_db_record("sg_snp_compare", {"_id": update_id}, {"vcf_path": filter_vcf_path})

    def add_sg_snp_compare_filter(self, vcf_path, task_id, project_sn, name=None):
        """
        新变异位点数据表
        :param vcf_path:
        :param task_id:
        :param project_sn:
        :return:
        """
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "vcf_path": vcf_path,
            "name": name if name else "pop_final_vcf",
            "desc": "",
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_snp_compare_filter", data_list)
        self.update_db_record("sg_snp_compare_filter", {"_id": main_id}, {"main_id": main_id})
        return main_id

if __name__ == "__main__":
    data = SnpCompare(None)
    task_id = "test"
    params = "无参导表"
    # vcf_path = "/mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/GeneticEvolution/pop.filtered.vcf"
    # data.add_sg_snp_compare(task_id, params)
    # data.add_sg_snp_compare_stat("/mnt/ilustre/users/sanger-dev/sg-users/zhaobinbin/wgs_noref/api/snp_compare/snp_compare", "5c1dcc44a4e1af553dd0f24b")
    data.add_snp_compare_detail("/mnt/ilustre/users/sanger-dev/workspace/20190125/SnpCompare_tsg_33313_0125173254851956_6936/output/snp_compare/snp_compare.table.xls", "5c1dcc44a4e1af553dd0f24b", "lala")
    # data.add_sg_snp_compare_filter(vcf_path, "test", "")



