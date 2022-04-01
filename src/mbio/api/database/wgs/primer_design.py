# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.23

from api_base import ApiBase
import datetime
import os


class PrimerDesign(ApiBase):
    def __init__(self, bind_object):
        """
        WGS 引物分析导表
        """
        super(PrimerDesign, self).__init__(bind_object)
        self._project_type = "dna_wgs"

    def add_sg_primer(self, project_sn, task_id, params=None, name=None):
        """
        sg_primer
        """
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_prinmer_design",
            "params": params if params else "null",
            "desc": "引物设计分析主表"
        }
        data_list.append(insert_data)
        main_id = self.db["sg_primer"].insert_one(insert_data).inserted_id
        self.update_db_record("sg_primer", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_primer_detail(self, primer_id, primer_result, download_file=None):
        """
        sg_primer_detail
        primer_result: variation.result
        """
        primer_id = self.check_objectid(primer_id)
        self.check_exists(primer_result)
        if download_file:
            self.update_db_record("sg_primer", {"_id": primer_id}, {"download_path": download_file})
        data_list = []
        num = 0
        with open(primer_result, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                try:
                    insert_data = {
                        "primer_id": primer_id,
                        "chr": item[0],
                        "pos": int(item[1]),
                        "ref": item[2],
                        "alt": item[3],
                        "forward_primer": item[4],
                        "forward_tm": item[5],
                        "forward_gc": item[6],
                        "forward_len": item[7],
                        "reverse_primer": item[8],
                        "reverse_tm": item[9],
                        "reverse_gc": item[10],
                        "reverse_len": item[11],
                        "product_size": item[12],
                        "variation_start": item[13],
                        "variation_end": item[14]
                    }
                except:
                    insert_data = {
                        "primer_id": primer_id,
                        "chr": item[0],
                        "pos": int(item[1]),
                        "ref": item[2],
                        "alt": item[3],
                        "forward_primer": "-",
                        "forward_tm": "-",
                        "forward_gc": "-",
                        "forward_len": "-",
                        "reverse_primer": "-",
                        "reverse_tm": "-",
                        "reverse_gc": "-",
                        "reverse_len": "-",
                        "product_size": "-",
                        "variation_start": "-",
                        "variation_end": "-"
                    }
                data_list.append(insert_data)
                # try:
                #     num = len(item[4].split(";"))
                #     for i in range(num):
                #         insert_data["forward_primer"] = item[4].split(";")[i]
                #         insert_data["forward_tm"] = float(item[5].split(";")[i])
                #         insert_data["forward_gc"] = float(item[6].split(";")[i])
                #         insert_data["forward_len"] = int(item[7].split(";")[i])
                #         insert_data["reverse_primer"] = item[8].split(";")[i]
                #         insert_data["reverse_tm"] = float(item[9].split(";")[i])
                #         insert_data["reverse_gc"] = float(item[10].split(";")[i])
                #         insert_data["reverse_len"] = int(item[11].split(";")[i])
                #         insert_data["product_size"] = int(item[12].split(";")[i])
                #         insert_data["variation_start"] = int(item[13].split(";")[i])
                #         insert_data["variation_end"] = int(item[14].split(";")[i])
                #         if "_id" in insert_data.keys():
                #             insert_data.pop("_id")
                #         # data_list.append(insert_data)
                #         self.db["sg_primer_detail"].insert_one(insert_data)
                # except:
                #     for i in range(num):
                #         insert_data["forward_primer"] = "-"
                #         insert_data["forward_tm"] = "-"
                #         insert_data["forward_gc"] = "-"
                #         insert_data["forward_len"] = "-"
                #         insert_data["reverse_primer"] = "-"
                #         insert_data["reverse_tm"] = "-"
                #         insert_data["reverse_gc"] = "-"
                #         insert_data["reverse_len"] = "-"
                #         insert_data["product_size"] = "-"
                #         insert_data["variation_start"] = "-"
                #         insert_data["variation_end"] = "-"
                #         if "_id" in insert_data.keys():
                #             insert_data.pop("_id")
                #         # data_list.append(insert_data)
                #         self.db["sg_primer_detail"].insert_one(insert_data)
        # self.db["sg_primer_detail"].insert_many(data_list)
        self.col_insert_data("sg_primer_detail", data_list)

    def add_sg_primer_detail_new(self, primer_id, primer_result):
        """
        sg_primer_detail
        primer_result: variation.result
        """
        primer_id = self.check_objectid(primer_id)
        self.check_exists(primer_result)
        data_list = []
        with open(primer_result, "r") as f:
            head = f.readline().strip().split("\t")
            num = len(head[4:]) / 11
            for line in f:
                item = line.strip().split("\t")
                forward_primer, forward_tm, forward_gc, forward_len = [], [], [], []
                reverse_primer, reverse_tm, reverse_gc, reverse_len = [], [], [], []
                product_size, variation_start, variation_end = [], [], []
                for i in range(num):
                    j = 11 * i
                    try:
                        forward_primer.append(item[4+j])
                        forward_tm.append(item[5+j])
                        forward_gc.append(item[6+j])
                        forward_len.append(item[7+j])
                        reverse_primer.append(item[8+j])
                        reverse_tm.append(item[9+j])
                        reverse_gc.append(item[10+j])
                        reverse_len.append(item[11+j])
                        product_size.append(item[12+j])
                        variation_start.append(item[13+j])
                        variation_end.append(item[14+j])
                    except:
                        break
                insert_data = {
                    "primer_id": primer_id,
                    "chr": item[0],
                    "pos": int(item[1]),
                    "ref": item[2],
                    "alt": item[3],
                    "forward_primer": ";".join(forward_primer),
                    "forward_tm": ";".join(forward_tm),
                    "forward_gc": ";".join(forward_gc),
                    "forward_len": ";".join(forward_len),
                    "reverse_primer": ";".join(reverse_primer),
                    "reverse_tm": ";".join(reverse_tm),
                    "reverse_gc": ";".join(reverse_gc),
                    "reverse_len": ";".join(reverse_len),
                    "product_size": ";".join(product_size),
                    "variation_start": ";".join(variation_start),
                    "variation_end": ";".join(variation_end)
                }
                data_list.append(insert_data)
        self.col_insert_data("sg_primer_detail", data_list)


if __name__ == "__main__":
    a = PrimerDesign(None)
    project_sn = 'wgs_test'
    task_id = 'wgs_test'
    # primer_id = a.add_sg_primer(project_sn, task_id)
    # primer_id = "5af53fc0a4e1af1416ba52e4"
    # primer_result = "/mnt/ilustre/users/sanger-dev/workspace/20180511/Single_wgs_test_0511150120_1273_6942/PrimerDesign/output/variation.result"
    # primer_result = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/test_file/primer/variation.result"
    # a.add_sg_primer_detail(primer_id, primer_result)
    primer_id = "5b1a04bf77b3f3b113e19de8"
    primer_result = "/mnt/ilustre/users/sanger-dev/workspace/20180608/Single_primer_design7/PrimerDesign/output/variant_result.xls"
    a.add_sg_primer_detail_new(primer_id, primer_result)
