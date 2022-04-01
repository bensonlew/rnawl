# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.05.18

from api_base import ApiBase
import datetime
import os


class SsrCompare(ApiBase):
    def __init__(self, bind_object):
        """
        WGS SSR分析导表
        """
        super(SsrCompare, self).__init__(bind_object)
        self._project_type = "dna_wgs"

    def add_sg_ssr_compare(self, project_sn, task_id, params=None, name=None):
        """
        sg_ssr_compare
        """
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_ssr_compare",
            "params": params if params else "null",
            "desc": "SSR比较分析主表"
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_ssr_compare", data_list)
        self.update_db_record("sg_ssr_compare", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_ssr_compare_stat(self, ssr_compare_id, ssr_stat):
        """
        sg_ssr_compare_stat
        ssr_stat: ssr.stat
        """
        ssr_compare_id = self.check_objectid(ssr_compare_id)
        self.check_exists(ssr_stat)
        data_list = []
        with open(ssr_stat, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "ssr_compare_id": ssr_compare_id,
                    "chr": item[0],
                    "ssr_num": int(item[1]),
                    "c": int(item[2]),
                    "c_star": int(item[3]),
                    "p1": int(item[4]),
                    "p2": int(item[5]),
                    "p3": int(item[6]),
                    "p4": int(item[7]),
                    "p5": int(item[8]),
                    "p6": int(item[9]),
                }
                data_list.append(insert_data)
        if data_list:
            self.col_insert_data("sg_ssr_compare_stat", data_list)
        else:
            self.bind_object.logger.info("ssr比较分析的stat结果为空")

    def add_sg_ssr_compare_detail(self, ssr_compare_id, ssr_detail):
        """
        sg_ssr_compare_detail
        ssr_detail: specimen.result
        """
        ssr_compare_id = self.check_objectid(ssr_compare_id)
        self.check_exists(ssr_detail)
        data_list = []
        num = 0
        with open(ssr_detail, "r") as f:
            lines = f.readlines()
            header = lines[0].strip().split("\t")
            # specimen_ids = {"sample1": header[7], "sample2": header[8]}
            sample1 = header[7]
            sample2 = header[8]
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "ssr_compare_id": ssr_compare_id,
                    "chr": item[0],
                    "ssr_nr": item[1],
                    "ssr_type": item[2],
                    "ssr": item[3],
                    "size": item[4],
                    "start": int(item[5]),
                    "end": int(item[6]),
                    sample1: item[7],
                    sample2: item[8]
                }
                origin_id = self.db["sg_ssr_compare_detail"].insert_one(insert_data).inserted_id
                if item[7] == "T":
                    try:
                        num = len(item[9].split(";"))
                        for i in range(num):
                            primer_data = {
                                "ssr_id": ssr_compare_id,
                                "origin_id": origin_id,
                                "specimen_id": sample1,
                                "chr": item[0],
                                "ssr_type": item[2],
                                "size": item[4],
                                "start": int(item[5]),
                                "end": int(item[6]),
                                "forward_primer": item[9].split(";")[i],
                                "forward_tm": float(item[10].split(";")[i]),
                                "forward_gc": float(item[11].split(";")[i]),
                                "forward_len": int(item[12].split(";")[i]),
                                "reverse_primer": item[13].split(";")[i],
                                "reverse_tm": float(item[14].split(";")[i]),
                                "reverse_gc": float(item[15].split(";")[i]),
                                "reverse_len": int(item[16].split(";")[i]),
                                "product_size": int(item[17].split(";")[i])
                            }
                            data_list.append(primer_data)
                    except:
                        primer_data = {
                            "ssr_id": ssr_compare_id,
                            "origin_id": origin_id,
                            "specimen_id": sample1,
                            "chr": item[0],
                            "ssr_type": item[2],
                            "size": item[4],
                            "start": int(item[5]),
                            "end": int(item[6]),
                            "forward_primer": "",
                            "forward_tm": "",
                            "forward_gc": "",
                            "forward_len": "",
                            "reverse_primer": "",
                            "reverse_tm": "",
                            "reverse_gc": "",
                            "reverse_len": "",
                            "product_size": ""
                        }
                        data_list.append(primer_data)
                else:
                    primer_data = {
                        "ssr_id": ssr_compare_id,
                        "origin_id": origin_id,
                        "specimen_id": sample1,
                        "chr": item[0],
                        "ssr_type": item[2],
                        "size": item[4],
                        "start": int(item[5]),
                        "end": int(item[6]),
                        "forward_primer": item[9],
                        "forward_tm": item[10],
                        "forward_gc": item[11],
                        "forward_len": item[12],
                        "reverse_primer": item[13],
                        "reverse_tm": item[14],
                        "reverse_gc": item[15],
                        "reverse_len": item[16],
                        "product_size": item[17]
                    }
                    data_list.append(primer_data)
                if item[8] == "T":
                    try:
                        num = len(item[18].split(";"))
                        for i in range(num):
                            primer_data = {
                                "ssr_id": ssr_compare_id,
                                "origin_id": origin_id,
                                "specimen_id": sample1,
                                "chr": item[0],
                                "ssr_type": item[2],
                                "size": item[4],
                                "start": int(item[5]),
                                "end": int(item[6]),
                                "forward_primer": item[18].split(";")[i],
                                "forward_tm": float(item[19].split(";")[i]),
                                "forward_gc": float(item[20].split(";")[i]),
                                "forward_len": int(item[21].split(";")[i]),
                                "reverse_primer": item[22].split(";")[i],
                                "reverse_tm": float(item[23].split(";")[i]),
                                "reverse_gc": float(item[24].split(";")[i]),
                                "reverse_len": int(item[25].split(";")[i]),
                                "product_size": int(item[26].split(";")[i])
                            }
                            data_list.append(primer_data)
                    except:
                        primer_data = {
                            "ssr_id": ssr_compare_id,
                            "origin_id": origin_id,
                            "specimen_id": sample1,
                            "chr": item[0],
                            "ssr_type": item[2],
                            "size": item[4],
                            "start": int(item[5]),
                            "end": int(item[6]),
                            "forward_primer": "",
                            "forward_tm": "",
                            "forward_gc": "",
                            "forward_len": "",
                            "reverse_primer": "",
                            "reverse_tm": "",
                            "reverse_gc": "",
                            "reverse_len": "",
                            "product_size": ""
                        }
                        data_list.append(primer_data)
                else:
                    primer_data = {
                        "ssr_id": ssr_compare_id,
                        "origin_id": origin_id,
                        "specimen_id": sample1,
                        "chr": item[0],
                        "ssr_type": item[2],
                        "size": item[4],
                        "start": int(item[5]),
                        "end": int(item[6]),
                        "forward_primer": item[18],
                        "forward_tm": item[19],
                        "forward_gc": item[20],
                        "forward_len": item[21],
                        "reverse_primer": item[22],
                        "reverse_tm": item[23],
                        "reverse_gc": item[24],
                        "reverse_len": item[25],
                        "product_size": item[26]
                    }
                    data_list.append(primer_data)
        if data_list:
            self.col_insert_data("sg_ssr_primer", data_list)
        else:
            self.bind_object.logger.info("ssr比较分析的detail结果为空")

    def add_sg_ssr_compare_detail_new(self, ssr_compare_id, ssr_detail, download_file):
        """
        sg_ssr_compare_detail
        ssr_detail: specimen.result
        """
        ssr_compare_id = self.check_objectid(ssr_compare_id)
        self.check_exists(ssr_detail)
        data_list = []
        num = 0
        self.update_db_record("sg_ssr_compare", {"main_id": ssr_compare_id}, {"download_path": download_file})
        with open(ssr_detail, "r") as f:
            lines = f.readlines()
            header = lines[0].strip().split("\t")
            sample1 = header[7]
            sample2 = header[8]
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "ssr_compare_id": ssr_compare_id,
                    "chr": item[0],
                    "ssr_nr": item[1],
                    "ssr_type": item[2],
                    "ssr": item[3],
                    "size": item[4],
                    "start": int(item[5]),
                    "end": int(item[6]),
                    sample1: item[7],
                    sample2: item[8]
                }
                if item[7] == "T":
                    try:
                        insert_data["forward_primer"] = item[9]
                        insert_data["forward_tm"] = item[10]
                        insert_data["forward_gc"] = item[11]
                        insert_data["forward_len"] = item[12]
                        insert_data["reverse_primer"] = item[13]
                        insert_data["reverse_tm"] = item[14]
                        insert_data["reverse_gc"] = item[15]
                        insert_data["reverse_len"] = item[16]
                        insert_data["product_size"] = item[17]
                    except:
                        insert_data["forward_primer"] = ""
                        insert_data["forward_tm"] = ""
                        insert_data["forward_gc"] = ""
                        insert_data["forward_len"] = ""
                        insert_data["reverse_primer"] = ""
                        insert_data["reverse_tm"] = ""
                        insert_data["reverse_gc"] = ""
                        insert_data["reverse_len"] = ""
                        insert_data["product_size"] = ""
                else:
                    insert_data["forward_primer"] = item[9]
                    insert_data["forward_tm"] = item[10]
                    insert_data["forward_gc"] = item[11]
                    insert_data["forward_len"] = item[12]
                    insert_data["reverse_primer"] = item[13]
                    insert_data["reverse_tm"] = item[14]
                    insert_data["reverse_gc"] = item[15]
                    insert_data["reverse_len"] = item[16]
                    insert_data["product_size"] = item[17]
                if item[8] == "T":
                    try:
                        insert_data["forward_primer"] = item[18]
                        insert_data["forward_tm"] = item[19]
                        insert_data["forward_gc"] = item[20]
                        insert_data["forward_len"] = item[21]
                        insert_data["reverse_primer"] = item[22]
                        insert_data["reverse_tm"] = item[23]
                        insert_data["reverse_gc"] = item[24]
                        insert_data["reverse_len"] = item[25]
                        insert_data["product_size"] = item[26]
                    except:
                        insert_data["forward_primer"] = ""
                        insert_data["forward_tm"] = ""
                        insert_data["forward_gc"] = ""
                        insert_data["forward_len"] = ""
                        insert_data["reverse_primer"] = ""
                        insert_data["reverse_tm"] = ""
                        insert_data["reverse_gc"] = ""
                        insert_data["reverse_len"] = ""
                        insert_data["product_size"] = ""
                else:
                    insert_data["forward_primer"] = item[18]
                    insert_data["forward_tm"] = item[19]
                    insert_data["forward_gc"] = item[20]
                    insert_data["forward_len"] = item[21]
                    insert_data["reverse_primer"] = item[22]
                    insert_data["reverse_tm"] = item[23]
                    insert_data["reverse_gc"] = item[24]
                    insert_data["reverse_len"] = item[25]
                    insert_data["product_size"] = item[26]
                data_list.append(insert_data)
        if data_list:
            self.col_insert_data("sg_ssr_compare_detail", data_list)
        else:
            self.bind_object.logger.info("ssr比较分析的detail结果为空")


if __name__ == "__main__":
    a = SsrCompare(None)
    project_sn = "wgs_test"
    task_id = "wgs_test"
    # ssr_stat = "/mnt/ilustre/users/sanger-dev/workspace/20180521/Single_ssr_compare_stat1/SsrCompareStat/output/ssr_stat.xls"
    # ssr_detail = "/mnt/ilustre/users/sanger-dev/workspace/20180521/Single_ssr_compare_stat1/SsrCompareStat/output/ssr_detail.xls"
    # ssr_compare_id = a.add_sg_ssr_compare(project_sn, task_id, params=None, name=None)
    # a.add_sg_ssr_compare_stat(ssr_compare_id, ssr_stat)
    # a.add_sg_ssr_compare_detail(ssr_compare_id, ssr_detail)
    ssr_compare_id = "5b1b5cc977b3f3b11311da4b"
    ssr_stat = "/mnt/ilustre/users/sanger-test/workspace/20180608/SsrCompare_tsanger_30180_0608152431447894_4074/output/ssr_stat.xls"
    ssr_detail = "/mnt/ilustre/users/sanger-test/workspace/20180608/SsrCompare_tsanger_30180_0608152431447894_4074/output/ssr_detail.xls"
    download_file = "rerewrweset/files/m_188/188_5b03d16580da8/tsanger_30180/interaction_results/ssr_compare"
    a.add_sg_ssr_compare_stat(ssr_compare_id, ssr_stat)
    a.add_sg_ssr_compare_detail_new(ssr_compare_id, ssr_detail, download_file)
