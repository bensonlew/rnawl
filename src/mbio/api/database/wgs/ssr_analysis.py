# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.23

from api_base import ApiBase
import datetime
import os


class SsrAnalysis(ApiBase):
    def __init__(self, bind_object):
        """
        WGS SSR分析导表
        """
        super(SsrAnalysis, self).__init__(bind_object)
        self._project_type = "dna_wgs"

    def add_sg_ssr(self, project_sn, task_id, params=None, name=None):
        """
        sg_ssr
        """
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_ssr_ref",
            "params": params if params else "null",
            "desc": "SSR参考基因组主表"
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_ssr", data_list)
        self.update_db_record("sg_ssr", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_ssr_stat(self, ssr_id, ssr_stat):
        """
        sg_ssr_stat
        ssr_stat: ssr.stat
        """
        ssr_id = self.check_objectid(ssr_id)
        self.check_exists(ssr_stat)
        data_list = []
        with open(ssr_stat, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "ssr_id": ssr_id,
                    "chr": item[0],
                    "ssr_num": int(item[1]) if item[1] else '--',
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
        self.col_insert_data("sg_ssr_stat", data_list)

    def add_sg_ssr_detail(self, ssr_id, ssr_detail):
        """
        sg_ssr_detail
        ssr_detail: ref.result
        """
        ssr_id = self.check_objectid(ssr_id)
        self.check_exists(ssr_detail)
        data_list = []
        num = 0
        with open(ssr_detail, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "ssr_id": ssr_id,
                    "chr": item[0],
                    "ssr_nr": item[1],
                    "ssr_type": item[2],
                    "ssr": item[3],
                    "size": item[4],
                    "start": int(item[5]),
                    "end": int(item[6])
                }
                origin_id = self.db["sg_ssr_detail"].insert_one(insert_data).inserted_id
                try:
                    num = len(item[7].split(";"))
                    for i in range(num):
                        primer_data = {
                            "ssr_id": ssr_id,
                            "origin_id": origin_id,
                            "chr": item[0],
                            "ssr_type": item[2],
                            "size": item[4],
                            "start": int(item[5]),
                            "end": int(item[6]),
                            "forward_primer": item[7].split(";")[i],
                            "forward_tm": float(item[8].split(";")[i]),
                            "forward_gc": float(item[9].split(";")[i]),
                            "forward_len": int(item[10].split(";")[i]),
                            "reverse_primer": item[11].split(";")[i],
                            "reverse_tm": float(item[12].split(";")[i]),
                            "reverse_gc": float(item[13].split(";")[i]),
                            "reverse_len": int(item[14].split(";")[i]),
                            "product_size": int(item[15].split(";")[i])
                        }
                        data_list.append(primer_data)
                        # self.db["sg_ssr_primer"].insert_one(primer_data)
                except:
                    for i in range(num):
                        primer_data = {
                            "ssr_id": ssr_id,
                            "origin_id": origin_id,
                            "chr": item[0],
                            "ssr_type": item[2],
                            "size": item[4],
                            "start": int(item[5]),
                            "end": int(item[6]),
                            "forward_primer": "-",
                            "forward_tm": "-",
                            "forward_gc": "-",
                            "forward_len": "-",
                            "reverse_primer": "-",
                            "reverse_tm": "-",
                            "reverse_gc": "-",
                            "reverse_len": "-",
                            "product_size": "-"
                        }
                        data_list.append(primer_data)
                        # self.db["sg_ssr_primer"].insert_one(primer_data)
        self.col_insert_data("sg_ssr_primer", data_list)

    def add_sg_ssr_detail_new(self, ssr_id, ssr_detail, download_file):
        """
        sg_ssr_detail
        ssr_detail: ref.result
        """
        ssr_id = self.check_objectid(ssr_id)
        self.check_exists(ssr_detail)
        data_list = []
        num = 0
        self.update_db_record("sg_ssr", {"main_id": ssr_id}, {"download_path": download_file})
        with open(ssr_detail, "r") as f:
            head = f.readline().strip().split("\t")
            num = len(head[6:]) / 9
            for line in f:
                item = line.strip().split("\t")
                forward_primer, forward_tm, forward_gc, forward_len = [], [], [], []
                reverse_primer, reverse_tm, reverse_gc, reverse_len = [], [], [], []
                product_size = []
                for i in range(num):
                    j = 9 * i
                    try:
                        forward_primer.append(item[7+j])
                        forward_tm.append(item[8+j])
                        forward_gc.append(item[9+j])
                        forward_len.append(item[10+j])
                        reverse_primer.append(item[11+j])
                        reverse_tm.append(item[12+j])
                        reverse_gc.append(item[13+j])
                        reverse_len.append(item[14+j])
                        product_size.append(item[15+j])
                    except:
                        break
                insert_data = {
                    "ssr_id": ssr_id,
                    "chr": item[0],
                    "ssr_nr": item[1],
                    "ssr_type": item[2],
                    "ssr": item[3],
                    "size": item[4],
                    "start": int(item[5]),
                    "end": int(item[6]),
                    "forward_primer": ";".join(forward_primer),
                    "forward_tm": ";".join(forward_tm),
                    "forward_gc": ";".join(forward_gc),
                    "forward_len": ";".join(forward_len),
                    "reverse_primer": ";".join(reverse_primer),
                    "reverse_tm": ";".join(reverse_tm),
                    "reverse_gc": ";".join(reverse_gc),
                    "reverse_len": ";".join(reverse_len),
                    "product_size": ";".join(product_size)
                }
                data_list.append(insert_data)
        self.col_insert_data("sg_ssr_detail", data_list)


if __name__ == "__main__":
    a = SsrAnalysis(None)
    project_sn = '188_5b03d16580da8'
    task_id = 'tsanger_30180'
    ssr_id = a.add_sg_ssr(project_sn, task_id)
    ssr_stat = "/mnt/ilustre/users/sanger-dev/app/database/dna_wgs_geneome/Arabidopsis_thaliana/NCBI/GCF_000001735.3/2011.05.11/ssr.stat"
    a.add_sg_ssr_stat(ssr_id, ssr_stat)
    ssr_detail = "/mnt/ilustre/users/sanger-dev/app/database/dna_wgs_geneome/Arabidopsis_thaliana/NCBI/GCF_000001735.3/2011.05.11/ssr.ref.result1.xls"
    download_file = ""
    a.add_sg_ssr_detail_new(ssr_id, ssr_detail, download_file)
