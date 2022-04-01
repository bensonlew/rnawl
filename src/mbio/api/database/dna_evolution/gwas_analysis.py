# -*- coding: utf-8 -*-
# __author__ = 'qingmei'
# modified 2018.0825

from api_base import ApiBase
from collections import defaultdict
import os
import datetime
import json
import math


class GwasAnalysis(ApiBase):
    """
    群体进化，GWAS关联分析接口
    """
    def __init__(self, bind_object):
        super(GwasAnalysis, self).__init__(bind_object)
        self._project_type = "dna_evolution"

    def add_sg_gwas_analysis(self, project_sn, task_id, upload_trit_path, params=None, name=None):
        """
        GWAS主表
        upload_trit_path:rere路径恢复成绝对路径
        module需要更新主表，增加字段csv_dir
        """
        if params:
            new_params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        else:
            new_params = "null"
        self.check_exists(upload_trit_path)
        trait_list = []
        with open(upload_trit_path, 'r') as r:
            line = r.readline()
            header = line.strip().split('\t')
            for trait in header[1:]:
                if trait in trait_list:
                    self.bind_object.logger.info("上传文件trait_id出现重复！{}".format(trait))
                trait_list.append(trait)
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "params": new_params,
            "name": name if name else "origin_sg_gwas_analysis",
            "desc": "GWAS关联分析接口",
            "trait_list": trait_list
        }
        main_id = self.db['sg_gwas_analysis'].insert_one(insert_data).inserted_id
        self.update_db_record("sg_gwas_analysis", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_gwas_analysis_stat(self, gwas_id, trait_stat, task_id, bar_dir):
        """
        性状细节表
        """
        gwas_id = self.check_objectid(gwas_id)   # 检查id是否是OBjectID
        self.check_exists(trait_stat)   # 检查文件是否存在
        # data_list = []
        with open(trait_stat, 'r') as r:
            lines = r.readlines()[1:]
            for line in lines:
                if line.startswith("#") or line.startswith("@"):
                    pass
                else:
                    tmp = line.strip().split('\t')
                    insert_data = {     # 文件内顺序
                        "gwas_id": gwas_id,
                        "trait_name": tmp[0],
                        "phenotype_num": tmp[1],
                        "average_value": tmp[2],
                        "variance": tmp[3],
                        "p_value": tmp[4],
                    }
                    # data_list.append(insert_data)
                    analysis_id = self.db["sg_gwas_analysis_stat"].insert_one(insert_data).inserted_id
                    self.add_sg_feature_bar(task_id, analysis_id, bar_dir, str(tmp[0]))
        # if len(data_list) == 0:  # add by hongdong 20180409
        #     self.bind_object.logger.info("trait文件为空！{}".format(trait_stat))
        # else:
        #     self.col_insert_data("sg_gwas_analysis_stat", data_list)
        print("sg_gwas_analysis_stat导入")

    def add_sg_feature_bar(self, task_id, gwas_id, bar_dir, trait_name):
        gwas_id = self.check_objectid(gwas_id)
        self.check_exists(bar_dir)
        for f in os.listdir(bar_dir):
            trit = f.split(".bar.txt")[0]
            if trit == trait_name:
                f_ = os.path.join(bar_dir, f)
                categories, value = [], []
                with open(f_, "r") as f1:
                    lines = f1.readlines()
                    for line in lines:
                        item = line.strip().split("\t")
                        categories.append(item[0])
                        value.append(int(item[1]))
                bar_id = self.add_sg_bar(gwas_id, task_id, categories)
                self.add_sg_bar_detail(bar_id, trait_name, value)
        print("sg_feature_bar导入")

    def add_sg_bar(self, feature_id, task_id, categories):
        """
        sg_bar
        """
        insert_data = {
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "name": "",
            "task_id": task_id,
            "origin_id": feature_id,
            "location": "gwas_feature_bar",
            "other_attr": "",
            "type": 4,
            "ext_name": "",
            "categories": categories
        }
        bar_id = self.db["sg_bar"].insert_one(insert_data).inserted_id
        print("sg_bar导入")
        return bar_id

    def add_sg_bar_detail(self, bar_id, name, values):
        """
        name:性状名称
        sg_bar_detail
        """
        bar_id = self.check_objectid(bar_id)
        insert_data = {
            "bar_id": bar_id,
            "name": name,
            "value": values
        }
        self.db["sg_bar_detail"].insert_one(insert_data)
        print("sg_bar_detail导入")

    def add_sg_manhattan(self, data_path, main_id, task_id, ext_name, chr_list, pos_list, value_list, trait, pvalue_default):
        """
        sg_manhattan
        data_path:/mnt/ilustre/users/sanger-dev/workspace/20181009/Single_gwas_analysis_12_20181009/GwasAnalysis/output/gwas_mvp/
        """
        origin_id = self.check_objectid(main_id)
        self.check_exists(data_path)
        if ext_name == "farmcpu":
            insert_data_1 = {
                "origin_id": origin_id,
                "task_id": task_id,
                "type": 2,
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "chr_list": chr_list,
                "name": str(trait),
                "pos": pos_list,
                "ext_name": "farmcpu",
                "location": "gwas_manhattan",
                "value": value_list,
                "pvalue_default": pvalue_default
            }
            manhattan_id_1 = self.db["sg_manhattan"].insert_one(insert_data_1).inserted_id
        if ext_name == "glm":
            insert_data_2 = {
                "origin_id": origin_id,
                "task_id": task_id,
                "type": 2,
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "chr_list": chr_list,
                "name": str(trait),
                "pos": pos_list,
                "ext_name": "glm",
                "location": "gwas_manhattan",
                "value": value_list,
                "pvalue_default": pvalue_default
            }
            manhattan_id_2 = self.db["sg_manhattan"].insert_one(insert_data_2).inserted_id
        if ext_name == "mlm":
            insert_data_3 = {
                "origin_id": origin_id,
                "task_id": task_id,
                "type": 2,
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "chr_list": chr_list,
                "name": str(trait),
                "pos": pos_list,
                "ext_name": "mlm",
                "location": "gwas_manhattan",
                "value": value_list,
                "pvalue_default": pvalue_default
            }
            manhattan_id_3 = self.db["sg_manhattan"].insert_one(insert_data_3).inserted_id

    # def add_sg_manhattan_detail(self, manhattan_id, name, x_categories, value):
    #     """
    #     sg_manhattan_detail
    #     """
    #     manhattan_id = self.check_objectid(manhattan_id)
    #     insert_data = {
    #         "manhattan_id": manhattan_id,
    #         "name": name,
    #         "pos_data": x_categories,
    #         "value": value
    #     }
    #     self.db["sg_manhattan_detail"].insert_one(insert_data)

    def add_sg_scatter(self, task_id, origin_id, gwas_output_dir, x_max, y_max, farmcpu_value, mlm_value, glm_value, trait):
        """
        sg_scatter
        """
        origin_id = self.check_objectid(origin_id)
        self.check_exists(gwas_output_dir)
        insert_data = {
            "origin_id": origin_id,
            "task_id": task_id,
            "name": "Gwas期望分布图",
            "type": 3,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "x_max": x_max,
            "y_max": y_max,
            "trait": trait,
            "location": "gwas_scatter"
        }
        scatter_id = self.db["sg_scatter"].insert_one(insert_data).inserted_id
        self.add_sg_scatter_detail(scatter_id, farmcpu_value, "FarmCPU")
        self.add_sg_scatter_detail(scatter_id, mlm_value, "MLM")
        self.add_sg_scatter_detail(scatter_id, glm_value, "GLM")

    def add_sg_scatter_detail(self, scatter_id, value, name):
        """
        sg_scatter_detail
        """
        scatter_id = self.check_objectid(scatter_id)
        insert_data = {
            "scatter_id": scatter_id,
            "name": name,
            "value": value
        }
        self.db["sg_scatter_detail"].insert_one(insert_data)

    def update_sg_gwas_analysis(self, gwas_id, csv_dir):
        """
        更新主表字段csv_dir用
        """
        self.update_db_record("sg_gwas_analysis", {"_id": gwas_id}, {"csv_dir": csv_dir})

    def check_exists(self, file_path):
        """
        用于检查文件及文件夹是否存在
        :param file_path:
        :return:
        """
        if not os.path.exists(file_path):
            raise Exception("文件或文件夹{}不存在！".format(file_path))

    def splist(self, l, s):
        """
        用于自动等分list
        :param l:list
        :param s:按多少等分
        :return:
        """
        return [l[i:i + s] for i in range(len(l)) if i % s == 0]


if __name__ == "__main__":
    a = GwasAnalysis(None)
    # project_sn = 'dna_evolution_cui'
    task_id = 'test_liu'
    # params = "{\"\":\"\"}"
    # # trait_list = ["t1", "t2", "t3", "t4"]
    # upload_trit_path = "/mnt/ilustre/users/sanger-dev/home/cuiqingmei/1.project/4.evolution/2.TraitAnnovar/trit3.xls"
    # gwas_id = a.add_sg_gwas_analysis(project_sn, task_id, upload_trit_path, params)
    # # summary_path接口task_id查表(tool不支持查表)，传给tool
    # trit_stat = "/mnt/ilustre/users/sanger-dev/workspace/20180822/Single_region_anno_0817_220180822/TraitAnnovar/output/annovar.trait.xls"
    # a.add_sg_gwas_analysis_stat(gwas_id, trit_stat)
    # bar_dir = "/mnt/ilustre/users/sanger-dev/workspace/20180822/Single_region_anno_0817_220180822/TraitAnnovar/bar_dir"
    # a.add_sg_feature_bar(task_id, gwas_id, bar_dir)
    origin_id = "5b836741a4e1af6072e8712d"
    gwas_id = "5b836741a4e1af6072e8712d"
    trait_stat = "/mnt/ilustre/users/sanger-dev/workspace/20181009/Single_gwas_analysis_12_20181009/GwasAnalysis/output/trait_annovar/annovar.trait.xls"
    task_id = "test123456"
    bar_dir = "/mnt/ilustre/users/sanger-dev/workspace/20181009/Single_gwas_analysis_12_20181009/GwasAnalysis/TraitAnnovar/bar_dir"
    a.add_sg_gwas_analysis_stat(gwas_id, trait_stat, task_id, bar_dir)
    # data_path = "/mnt/ilustre/users/sanger-dev/workspace/20181009/Single_gwas_analysis_12_20181009/GwasAnalysis/output/gwas_mvp/"
    # gwas_output_dir = "/mnt/ilustre/users/sanger-dev/workspace/20181015/Single_gwas_analysis_7_20181015/GwasAnalysis/output/gwas_mvp"
    # a.add_sg_manhattan(data_path, main_id)
    # a.add_sg_scatter(task_id, origin_id, gwas_output_dir)
"""
已完成导表书写：
    sg_gwas_analysis：GWAS关联分析接口
    sg_gwas_analysis_stat：性状分布统计表
    sg_bar、sg_bar_detail：性状分布统计表--关联的直方图主表和细节表
"""
"""
未导表：
    sg_?：期望分布图主表
    sg_?_detail：期望分布图细节表
"""
