# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 20181220

import os
import math
import types
import datetime
from bson.son import SON
from bson.objectid import ObjectId
from biocluster.api.database.base import Base, report_check


class SpecimenQc(Base):
    """
    无参WGS导表：软件列表、样本信息、质控
    """
    def __init__(self, bind_object):
        super(SpecimenQc, self).__init__(bind_object)
        self._project_type = "dna_noref_wgs"

    def check_exists(self, file):
        """
        检查file文件是否存在
        """
        if not os.path.exists(file):
            raise Exception("文件：%s不存在，请检查" % file)
        return True

    def check_objectid(self, id):
        """
        检查id是否是mongo的_id
        """
        if not isinstance(id, ObjectId):
            if isinstance(id, types.StringTypes):
                id = ObjectId(id)
            else:
                raise Exception("id:%s必须为ObjectID或其对应的字符串" % id)
        return id

    def add_sg_task(self, member_id, member_type, cmd_id, project_sn, task_id):
        """
        sg_task
        """
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "member_id": member_id,
            "member_type": member_type,
            "cmd_id": cmd_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        main_id = self.db["sg_task"].insert_one(insert_data).inserted_id
        print "add_sg_task success!"

    def add_sg_software(self, task_id, item, software_name, version, params=None):
        """
        sg_software
        """
        insert_data = {
            "task_id": task_id,
            "item": item,
            "software_name": software_name,
            "version": version,
            "params": params,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        self.db["sg_software"].insert_one(insert_data)
        print "sg_software success!"

    def add_sg_specimen(self, task_id, initial_name, analysis_name, library, raw_data, clean_data, desc=None):
        """
        sg_specimen
        """
        insert_data = {
            "task_id": task_id,
            "initial_name": initial_name,
            "analysis_name": analysis_name,
            "new_name": None,
            "library": library,
            "batch": "batch1",
            "raw_data": raw_data,
            "clean_data": clean_data,
            "desc": desc if desc else "sg_specimen",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        specimen_id = self.db["sg_specimen"].insert_one(insert_data).inserted_id
        print "add_sg_specimen success!"
        return specimen_id

    def add_sg_specimen_qc(self, task_id, qc_stat, raw_dir, clean_dir):
        """
        sg_specimen_qc
        """
        self.check_exists(qc_stat)
        with open(qc_stat, "r") as f:
            for line in f:
                if line.startswith("#"):
                    pass
                else:
                    item = line.strip().split("\t")
                    if item[0] in ["1080", "1145", "1828"]:
                        specimen_id = "a" + item[0]
                    else:
                        specimen_id = item[0]
                    insert_data = {
                        "task_id": task_id,
                        "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        "specimen_id": specimen_id,
                        "raw_reads": int(item[1]),
                        "raw_base": int(item[2]),
                        "raw_gc_rate": float(item[3]),
                        "raw_q30_rate": float(item[4]),
                        "clean_reads": int(item[5]),
                        "clean_base": int(item[6]),
                        "clean_gc_rate": float(item[7]),
                        "clean_q30_rate": float(item[8])
                    }
                    main_id = self.db["sg_specimen_qc"].insert_one(insert_data).inserted_id
                    self.db["sg_specimen_qc"].update_one({"_id": main_id}, {"$set": {"main_id": main_id}})
                    raw_atgc_path = os.path.join(raw_dir, "atgc/" + item[0] + ".atgc.xls")
                    clean_atgc_path = os.path.join(clean_dir, "atgc/" + item[0] + ".atgc.xls")
                    raw_qual_path = os.path.join(raw_dir, "qual/" + item[0] + ".qual.xls")
                    clean_qual_path = os.path.join(clean_dir, "qual/" + item[0] + ".qual.xls")
                    a.add_sg_base_content(task_id=task_id, specimen_qc_id=main_id, qual_path=raw_atgc_path, location="raw_reads", sample=specimen_id)
                    a.add_sg_base_content(task_id=task_id, specimen_qc_id=main_id, qual_path=clean_atgc_path, location="clean_reads", sample=specimen_id)
                    a.add_sg_base_error(task_id=task_id, specimen_qc_id=main_id, qual_path=raw_qual_path, location="error_raw_reads", sample=specimen_id)
                    a.add_sg_base_error(task_id=task_id, specimen_qc_id=main_id, qual_path=clean_qual_path, location="error_clean_reads", sample=specimen_id)
        print "add_sg_specimen_qc success!"

    def add_sg_curve(self, task_id, specimen_qc_id, name, categories, location, sample, types=1):
        if not isinstance(specimen_qc_id, ObjectId):
            if isinstance(specimen_qc_id, types.StringTypes):
                specimen_qc_id = ObjectId(specimen_qc_id)
            else:
                raise Exception("specimen_qc_id必须为ObjectID或其对应的字符串")
        insert_data = {
            "task_id": task_id,
            "origin_id": specimen_qc_id,
            "name": name,
            "categories": categories,
            "type": types,
            "location": location,
            "other_attr": "",
            "ext_name": {"type": "sample", "title": sample},
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        coll = self.db["sg_curve"]
        main_id = coll.insert_one(insert_data).inserted_id
        print "add_sg_curve success!"
        return main_id

    def add_sg_curve_detail(self, curve_id, name, value):
        if not isinstance(curve_id, ObjectId):
            if isinstance(curve_id, types.StringTypes):
                curve_id = ObjectId(curve_id)
            else:
                raise Exception("curve_id必须为ObjectID或其对应的字符串")
        insert_data = {
            "curve_id": curve_id,
            "name": name,
            "value": value
        }
        coll = self.db["sg_curve_detail"]
        main_id = coll.insert_one(insert_data).inserted_id
        print "add_sg_curve_detail success!"

    def add_sg_base_content(self, task_id, specimen_qc_id, qual_path, location, sample):
        """
        碱基质量分布导表
        """
        if not isinstance(specimen_qc_id, ObjectId):
            if isinstance(specimen_qc_id, types.StringTypes):
                specimen_qc_id = ObjectId(specimen_qc_id)
            else:
                raise Exception("specimen_qc_id必须为ObjectID或其对应的字符串")
        if not os.path.exists(qual_path):
            raise Exception("{}文件路径不存在，请检查!".format(qual_path))
        categories = []
        a_list, t_list, g_list, c_list, n_list = [], [], [], [], []
        with open(qual_path, "r") as f:
            header = f.readline()
            for line in f:
                item = line.strip().split("\t")
                categories.append(item[0])
                a_list.append(None if math.isnan(float(item[1])) else float(item[1]))
                t_list.append(None if math.isnan(float(item[2])) else float(item[2]))
                g_list.append(None if math.isnan(float(item[3])) else float(item[3]))
                c_list.append(None if math.isnan(float(item[4])) else float(item[4]))
                n_list.append(None if math.isnan(float(item[5])) else float(item[5]))
        main_id = self.add_sg_curve(task_id=task_id, specimen_qc_id=specimen_qc_id, name=sample,\
                  categories=categories, location=location, sample=sample, types=1)
        self.add_sg_curve_detail(curve_id=main_id, name="A", value=a_list)
        self.add_sg_curve_detail(curve_id=main_id, name="T", value=t_list)
        self.add_sg_curve_detail(curve_id=main_id, name="C", value=c_list)
        self.add_sg_curve_detail(curve_id=main_id, name="G", value=g_list)
        self.add_sg_curve_detail(curve_id=main_id, name="N", value=n_list)

    def add_sg_base_error(self, task_id, specimen_qc_id, qual_path, location, sample):
        if not isinstance(specimen_qc_id, ObjectId):
            if isinstance(specimen_qc_id, types.StringTypes):
                specimen_qc_id = ObjectId(specimen_qc_id)
            else:
                raise Exception("specimen_qc_id必须为ObjectID或其对应的字符串")
        if not os.path.exists(qual_path):
            raise Exception("{}文件路径不存在，请检查!".format(qual_path))
        categories = []
        e_list = []
        # sample = os.path.basename(qual_path).split(".")[0]
        with open(qual_path, "r") as f:
            header = f.readline()
            for line in f:
                item = line.strip().split("\t")
                categories.append(item[0])
                if math.isnan(float(item[6])):
                    error = None
                else:
                    error = 1 / (10 ** (float(item[-1]) / 10))
                e_list.append(error)
        main_id = self.add_sg_curve(task_id=task_id, specimen_qc_id=specimen_qc_id, name="error",\
                  categories=categories, location=location, sample=sample, types=1)
        self.add_sg_curve_detail(curve_id=main_id, name="error", value=e_list)



if __name__ == "__main__":
    a = SpecimenQc(None)
    project_sn = "188_5c039c5db67e0"
    task_id = "noref_test1"
    # a.add_sg_task(member_id="m_188", member_type=1, cmd_id=196, project_sn=project_sn, task_id=task_id)
    # a.add_sg_software(task_id=task_id, item="SNP变异检测", software_name="ipyrad", version="0.01")
    qc_stat = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/qc.txt"
    raw_dir = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/01.fastq_qc/rawdata_qc"
    clean_dir = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/01.fastq_qc/cleandata_qc"
    a.add_sg_specimen_qc(task_id, qc_stat, raw_dir, clean_dir)
    # sample = "1080"
    # specimen_qc_id = a.add_sg_specimen(task_id=task_id, initial_name=sample, analysis_name=sample, library="HCF120005", raw_data="1080.raw.1.fastq.gz;1080.raw.2.fastq.gz", clean_data="1080.clean.1.fastq.gz;1080.clean.2.fastq.gz")
    # atgc_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/01.fastq_qc/rawdata_qc/atgc/1080.atgc.xls"
    # a.add_sg_base_content(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=atgc_path, location="raw_reads", sample=sample)
    # atgc_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/01.fastq_qc/cleandata_qc/atgc/1080.atgc.xls"
    # a.add_sg_base_content(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=atgc_path, location="clean_reads", sample=sample)
    # qual_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/01.fastq_qc/rawdata_qc/qual/1080.qual.xls"
    # a.add_sg_base_error(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=qual_path, location="error_raw_reads", sample=sample)
    # qual_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/01.fastq_qc/cleandata_qc/qual/1080.qual.xls"
    # a.add_sg_base_error(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=qual_path, location="error_clean_reads", sample=sample)
    # sample = "1145"
    # specimen_qc_id = a.add_sg_specimen(task_id=task_id, initial_name=sample, analysis_name=sample, library="HCF120006", raw_data="1145.raw.1.fastq.gz;1145.raw.2.fastq.gz", clean_data="1145.clean.1.fastq.gz;1145.clean.2.fastq.gz")
    # atgc_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/01.fastq_qc/rawdata_qc/atgc/1145.atgc.xls"
    # a.add_sg_base_content(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=atgc_path, location="raw_reads", sample=sample)
    # atgc_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/01.fastq_qc/cleandata_qc/atgc/1145.atgc.xls"
    # a.add_sg_base_content(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=atgc_path, location="clean_reads", sample=sample)
    # qual_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/01.fastq_qc/rawdata_qc/qual/1145.qual.xls"
    # a.add_sg_base_error(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=qual_path, location="error_raw_reads", sample=sample)
    # qual_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/01.fastq_qc/cleandata_qc/qual/1145.qual.xls"
    # a.add_sg_base_error(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=qual_path, location="error_clean_reads", sample=sample)
    # sample = "c29875"
    # specimen_qc_id = a.add_sg_specimen(task_id=task_id, initial_name=sample, analysis_name=sample, library="HCF120044", raw_data="c29875.raw.1.fastq.gz;c29875.raw.2.fastq.gz", clean_data="c29875.clean.1.fastq.gz;c29875.clean.2.fastq.gz")
    # atgc_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/01.fastq_qc/rawdata_qc/atgc/c29875.atgc.xls"
    # a.add_sg_base_content(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=atgc_path, location="raw_reads", sample=sample)
    # atgc_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/01.fastq_qc/cleandata_qc/atgc/c29875.atgc.xls"
    # a.add_sg_base_content(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=atgc_path, location="clean_reads", sample=sample)
    # qual_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/01.fastq_qc/rawdata_qc/qual/c29875.qual.xls"
    # a.add_sg_base_error(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=qual_path, location="error_raw_reads", sample=sample)
    # qual_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/01.fastq_qc/cleandata_qc/qual/c29875.qual.xls"
    # a.add_sg_base_error(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=qual_path, location="error_clean_reads", sample=sample)
    # sample = "1828"
    # specimen_qc_id = a.add_sg_specimen(task_id=task_id, initial_name=sample, analysis_name=sample, library="HCF120045", raw_data="1828.raw.1.fastq.gz;1828.raw.2.fastq.gz", clean_data="1828.clean.1.fastq.gz;1828.clean.2.fastq.gz")
    # atgc_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/01.fastq_qc/rawdata_qc/atgc/1828.atgc.xls"
    # a.add_sg_base_content(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=atgc_path, location="raw_reads", sample=sample)
    # atgc_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/01.fastq_qc/cleandata_qc/atgc/1828.atgc.xls"
    # a.add_sg_base_content(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=atgc_path, location="clean_reads", sample=sample)
    # qual_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/01.fastq_qc/rawdata_qc/qual/1828.qual.xls"
    # a.add_sg_base_error(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=qual_path, location="error_raw_reads", sample=sample)
    # qual_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/noref_WGS/mongo_data/01.fastq_qc/cleandata_qc/qual/1828.qual.xls"
    # a.add_sg_base_error(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=qual_path, location="error_clean_reads", sample=sample)
