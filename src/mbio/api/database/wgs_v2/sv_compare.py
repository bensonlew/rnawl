# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190225

import os
import datetime
from types import StringTypes
from bson.objectid import ObjectId
from biocluster.api.database.base import Base, report_check


class SvCompare(Base):
    """
    sv比较分析导表
    """
    def __init__(self, bind_object):
        super(SvCompare, self).__init__(bind_object)
        self._project_type = "dna_wgs_v2"

    def add_sg_sv_compare(self, project_sn, task_id, subname=None, type="single", params=None):
        """
        sg_sv_compare导表
        """
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "status": "end",
            "params": params if params else "null",
            "type": type,
            "name": "SvCompareAnalysis",
            "desc": "SV比较分析",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        if type == "single":
            insert_data["subname"] = subname
        main_id = self.db["sg_sv_compare"].insert_one(insert_data).inserted_id
        self.db["sg_sv_compare"].update({"_id": main_id}, {"$set": {"main_id": main_id}})
        return main_id

    def add_sg_sv_compare_stat(self, compare_id, summary_path):
        """
        sg_sv_compare_stat导表
        summary_path: sample1_VS_sample2.summary.xls
        """
        if not isinstance(compare_id, ObjectId):
            if isinstance(compare_id, StringTypes):
                compare_id = ObjectId(compare_id)
            else:
                raise Exception("compare_id:%s必须为ObjectId对象或字符串，请检查" % compare_id)
        if not os.path.exists(summary_path):
            raise Exception("文件：%s不存在，请检查" % summary_path)
        sub_name = os.path.basename(summary_path).split(".summary.xls")[0]
        data_list = []
        with open(summary_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "compare_id": compare_id,
                    "name": sub_name,
                    "chr": item[0],
                    "del": int(item[1]),
                    "ins": int(item[2]),
                    "dup": int(item[3]),
                    "inv": int(item[4]),
                    "bnd": int(item[5])
                }
                data_list.append(insert_data)
        if data_list:
            self.db["sg_sv_compare_stat"].insert_many(data_list)
        # self.bind_object.logger.info("sg_sv_compare_stat导表成功：%s" % summary_path)

    def add_sg_sv_compare_detail(self, compare_id, detail_path, type="single"):
        """
        sg_sv_compare_detail导表
        detail_path: sample1_VS_sample2.detail.xls
        """
        if not isinstance(compare_id, ObjectId):
            if isinstance(compare_id, StringTypes):
                compare_id = ObjectId(compare_id)
            else:
                raise Exception("compare_id:%s必须为ObjectId对象或字符串，请检查" % compare_id)
        if not os.path.exists(detail_path):
            raise Exception("文件：%s不存在，请检查" % detail_path)
        sub_name = os.path.basename(detail_path).split(".detail.xls")[0]
        data_list = []
        with open(detail_path, "r") as f:
            lines = f.readlines()
            header = lines[0].strip().split("\t")
            length = len(header)
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "compare_id": compare_id,
                    "name": sub_name,
                    "sv_id": item[0],
                    "chr": item[1],
                    "start": int(item[2]),
                    "end": int(item[3]),
                    "len": int(item[4]),
                    "region": item[1] + ":" + item[2] + "-" + item[3]
                }
                for i in range(5, length):
                    insert_data[header[i]] = item[i]
                data_list.append(insert_data)
        if data_list:
            self.db["sg_sv_compare_detail"].insert_many(data_list)
        if type == "multiple":
            header_, header_dict = [], {}
            for h in header[5:]:
                if h.endswith("_genotype"):
                    header_.append(h.split("_genotype")[0] + " Genotype")
                    header_dict[h] = h.split("_genotype")[0] + " Genotype"
                if h.endswith("_depth"):
                    header_.append(h.split("_depth")[0] + " Allele Depth")
                    header_dict[h] = h.split("_depth")[0] + " Allele Depth"
            self.db["sg_sv_compare"].update_one({"main_id": compare_id}, {"$set": {"header": header_, "header_dict": header_dict}})
        # self.bind_object.logger.info("sg_sv_compare_detail导表成功：%s" % detail_path)

    def update_vcf_path(self, compare_id, output_dir, s3_output_dir):
        """
        更新主表的vcf_path
        """
        if not isinstance(compare_id, ObjectId):
            if isinstance(compare_id, StringTypes):
                compare_id = ObjectId(compare_id)
            else:
                raise Exception("compare_id:%s必须为ObjectId对象或字符串，请检查" % compare_id)
        if not os.path.exists(output_dir):
            raise Exception("文件：%s不存在，请检查" % output_dir)
        result = self.db["sg_sv_compare"].find_one({"main_id": compare_id})
        if result["type"] == "single":
            subname, vcf_path = [], []
            for f in os.listdir(output_dir):
                if f.endswith("detail.xls"):
                    vcf_file = os.path.join(s3_output_dir, f)
                    sub_name = os.path.basename(f).split(".detail.xls")[0]
                    subname.append(sub_name)
                    vcf_path.append(vcf_file)
            self.db["sg_sv_compare"].update_one({"main_id": compare_id}, {"$set": {"subname": subname, "vcf_path": vcf_path}})
        else:
            for f in os.listdir(output_dir):
                if f.endswith("detail.xls"):
                    vcf_path = os.path.join(s3_output_dir, f)
                    self.db["sg_sv_compare"].update_one({"main_id": compare_id}, {"$set": {"vcf_path": vcf_path}})


if __name__ == "__main__":
    a = SvCompare(None)
    project_sn = "wgs_v2"
    task_id = "wgs_v2"
    subname = ["SRR5739119_VS_SRR5739120"]
    summary_path = "/mnt/ilustre/users/sanger-dev/workspace/20190227/Single_sv_compare/SvCompare/output/SRR5739119_VS_SRR5739120.summary.xls"
    detail_path = "/mnt/ilustre/users/sanger-dev/workspace/20190227/Single_sv_compare/SvCompare/output/SRR5739119_VS_SRR5739120.detail.xls"
    compare_id = a.add_sg_sv_compare(project_sn, task_id, subname)
    a.add_sg_sv_compare_stat(compare_id, summary_path)
    a.add_sg_sv_compare_detail(compare_id, detail_path)
