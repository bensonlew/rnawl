# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.23

from api_base import ApiBase
from bson.objectid import ObjectId
import datetime
import os
import re


class VennAnalysis(ApiBase):
    def __init__(self, bind_object):
        """
        WGS venn分析导表
        """
        super(VennAnalysis, self).__init__(bind_object)
        self._project_type = "dna_wgs"

    def add_sg_venn_analysis(self, project_sn, task_id, params=None, name=None):
        """
        sg_venn_analysis
        """
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_venn_analysis",
            "params": params if params else "null",
            "desc": "venn分析主表"
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_venn_analysis", data_list)
        self.update_db_record("sg_venn_analysis", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def update_table_name(self, analysis_id, update_dict):
        """
        更新sg_venn_analysis，用于详情表的下拉框
        """
        main_id = self.check_objectid(analysis_id)
        self.update_db_record("sg_venn_analysis", {"main_id": main_id}, update_dict)
        print "更新sg_venn_analysis的table_name成功"

    def add_sg_venn_analysis_stat(self, venn_id, venn_stat):
        """
        sg_venn_analysis_stat
        venn_stat:
        """
        venn_id = self.check_objectid(venn_id)
        self.check_exists(venn_stat)
        data_list = []
        with open(venn_stat, "r") as f:
            lines = f.readlines()
            for line in lines:
                item = line.strip().split("\t")
                insert_data = {
                    "analysis_id": venn_id,
                    "type": item[0],
                    "num": int(item[1])
                }
                data_list.append(insert_data)
        if data_list:
            self.col_insert_data("sg_venn_analysis_stat", data_list)
        else:
            self.bind_object.logger.info("结果文件{}为空".format(venn_stat))

    def add_sg_venn_analysis_detail(self, venn_id, name, venn_detail):
        """
        sg_venn_analysis_detail
        venn_detail:
        """
        venn_id = self.check_objectid(venn_id)
        self.check_exists(venn_detail)
        data_list = []
        title = {}
        types = {}
        with open(venn_detail, "r") as f:
            lines = f.readlines()
            header = lines[0].strip().split("\t")
            for i in range(5, len(header)):
                title["key"+str(i-4)] = header[i]
                types["key"+str(i-4)] = i
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "analysis_id": venn_id,
                    "name": name,
                    "chr": item[0],
                    "pos": int(item[1]),
                    "ref": item[2],
                    "alt": item[3],
                    "anno": item[4]
                }
                for t in types.keys():
                    try:
                        insert_data[t] = item[types[t]]
                    except:
                        print line
                data_list.append(insert_data)
        self.col_insert_data("sg_venn_analysis_detail", data_list)
        # self.update_db_record("sg_venn_analysis", {"main_id": venn_id}, {"title": title})
        print "更新sg_venn_analysis的title成功"
        return title

    def add_sg_venn(self, venn_id, name=None):
        """
        sg_venn 画venn图主表
        """
        origin_id = self.check_objectid(venn_id)
        data_list = []
        insert_data = {
            "origin_id": origin_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "name": name if name else "venn",
            "type": 1,
            "location": "venn_analysis"
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_venn", data_list)
        self.update_db_record("sg_venn", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_venn_detail(self, venn_id, set_ids, diff_list):
        """
        sg_venn_detail
        diff_list: 本次分析的所有差异结果的list
        """
        venn_id = self.check_objectid(venn_id)
        attr = []
        data_list = []
        set_dict = {}
        for set_id in set_ids:
            result = self.col_find_one("sg_site_set", query_dic={"set_id": set_id})
            name = result["name"]
            set_dict[name] = set_id
        with open(diff_list, "r") as f:
            for line in f:
                item = line.strip().split("\t")
                venn_list = []
                with open(item[1], "r") as fp:
                    lines = fp.readlines()
                    for line in lines[1:]:
                        item1 = line.strip().split("\t")
                        str = item1[0] + "_" + item1[1] + "_" + item1[2] + "_" + item1[3]
                        venn_list.append(str)
                venn_list = list(set(venn_list))
                attr.append(item[0])
                insert_data = {
                    "venn_id": venn_id,
                    "name": set_dict[item[0]],
                    "venn_list": venn_list
                }
                data_list.append(insert_data)
                self.db["sg_venn_detail"].insert_one(insert_data)
        # self.col_insert_data("sg_venn_detail", data_list)
        self.update_db_record("sg_venn", {"_id": venn_id}, {"attr": attr})

    def get_variant_result_list(self, set_ids, type, diff_list, rerewrwese_path=None):
        """
        根据set_ids得到对应的比较分析结果文件
        """
        with open(diff_list, "w") as w:
            for set_id in set_ids:
                result = self.col_find_one("sg_site_set", query_dic={"set_id": set_id})
                if not result:
                    raise Exception("没有找到sg_site_set表中set_id为{}的记录，请检查".format(set_id))
                name = result["name"]
                if type == "snp":
                    diff_result = self.col_find_one("sg_snp_compare", query_dic={"main_id": ObjectId(set_id)})
                else:
                    diff_result = self.col_find_one("sg_indel_compare", query_dic={"main_id": ObjectId(set_id)})
                try:
                    if "download_path" in diff_result.keys():
                        if rerewrwese_path and not re.match('.*://.*', diff_result["download_path"]):  # 新旧数据兼容
                            diff_variant = os.path.join(rerewrwese_path, diff_result["download_path"].lstrip("/"))
                        else:
                            diff_variant = diff_result["download_path"]
                    else:
                        diff_variant = diff_result["diff_variant"]
                except:
                    raise Exception("根据set_id: 没找到对应的diff_variant".format(set_id))
                w.write(name + "\t" + diff_variant + "\n")


if __name__ == "__main__":
    a = VennAnalysis(None)
    project_sn = 'wgs_test'
    task_id = 'wgs_test'
    # venn_id = a.add_sg_venn_analysis(project_sn, task_id)
    # venn_id = "5add3e1ea4e1af2e0001c9a0"
    # a.add_sg_venn(venn_id)
    # venn_detail = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/test_file/venn/ipop.result"
    # a.add_sg_venn_analysis_stat(venn_id, venn_detail)
    venn_id = "5add3e1ea4e1af2e0001c9a0"
    venn_detail = "/mnt/ilustre/users/sanger-dev/workspace/20180511/Single_wgs_test_0511090312_1715_4265/VennAnalysis/output/GC_bulk_JY102_diff_SNP_20180507_103953.result"
    a.add_sg_venn_analysis_detail(venn_id, venn_detail)
    # venn_id = "5add3ec1a4e1af2e986e0089"
    # venn_paths = ["/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/test_file/venn/group1_vs_group2_diff_20180416_135419.list1",
    #               "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/test_file/venn/g1_vs_g2_diff_20180418_132846.list1"]
    # a.add_sg_venn_detail(venn_id, venn_paths)
    # set_ids = ['5ad933aba4e1af1e32509b67', '5ad93377a4e1af1e32509b65']
    # type = "snp"
    # diff_list = "/mnt/ilustre/users/sanger-dev/workspace/20180427/Single_venn_ana1/VennAnalysis/test.list"
    # a.get_variant_result_list(set_ids, type, diff_list)
