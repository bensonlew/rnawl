# -*- coding: utf-8 -*-
# __author__ = 'wentian.liu'
# modified 2019.04.10

from api_base import ApiBase
import datetime
import os
import json
from bson.objectid import ObjectId


class Circos(ApiBase):
    """
    circos图导表
    """
    def __init__(self, bind_object):
        super(Circos, self).__init__(bind_object)
        self._project_type = "dna_wgs_v2"

    def add_sg_circos(self, main_id, png, svg):
        """
        更新circos图主表
        """
        id = self.check_objectid(main_id)
        # self.check_exists(png)
        # self.check_exists(svg)
        task_id = self.col_find_one(collection="sg_circos", query_dic={"_id": id})['task_id']
        chr_list = self.col_find_one(collection="sg_circos",
                                     query_dic={"task_id": task_id, "name": "origin_circos"})['chr_list']
        self.update_db_record("sg_circos", {"_id": id}, {"circos_png_path": png})
        self.update_db_record("sg_circos", {"_id": id}, {"circos_svg_path": svg})
        self.update_db_record("sg_circos", {"_id": id}, {"chr_list": chr_list})

    def add_sg_origin_circos(self, task_id, task_type, submit_location, chromosome, color, variant, project_sn,
                             png, svg, total_chrlist):
        """
        添加一条name为origin_circos的记录用来给前端取出所有chr
        """
        chr_list = []
        with open(total_chrlist, "r")as fr:
            lines = fr.readlines()
            for line in lines:
                tmp = line.strip().split("\t")
                chr_list.append(tmp[0])
        params_json = {
            "task_type": int(task_type),
            "submit_location": submit_location,
            "task_id": task_id,
            "chromosome": chromosome,
            "color": color,
            "variant": variant,
            "color_type": "circle",
            "chongmingming_result": ""
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        insert_data = {
            "name": "origin_circos",
            "status": 'end',
            "project_sn": project_sn,
            "task_id": task_id,
            "params": params,
            "desc": "Circos主表！",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "circos_png_path": png,
            "circos_svg_path": svg,
            "chr_list": chr_list
        }
        main_id = self.db["sg_circos"].insert_one(insert_data).inserted_id
        self.update_db_record("sg_circos", {"_id": main_id}, {"main_id": main_id})


if __name__ == "__main__":
    a = Circos(None)
    task_id = "wgs_v2"
    task_type = 1
    submit_location = "test_wgs"
    chromosome = "chr1,chr2,chr3"
    color = 1
    variant = [{"variant":"snpplusindel","type":"before","analysis_object":"AH03","window":10000,"style":"histogram"}]
    status = "end"
    project_sn = "wgs_v2"
    png = "/mnt/ilustre/users/sanger-dev/workspace/20190410/Single_wgs_v2_0410133919871869_2239/Circos/output/circos.png"
    svg = "/mnt/ilustre/users/sanger-dev/workspace/20190410/Single_wgs_v2_0410133919871869_2239/Circos/output/circos.svg"
    total_chrlist = "/mnt/ilustre/users/sanger-dev/app/database/dna_wgs_geneome/Citrus_sinensis/HZAU/V2/2013.6.21/total.chrlist"
    a.add_sg_origin_circos(task_id, task_type, submit_location, chromosome, color, variant, status, project_sn, png, svg, total_chrlist)