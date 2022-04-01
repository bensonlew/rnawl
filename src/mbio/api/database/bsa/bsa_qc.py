# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.02.06


# import pymongo
import os
import types
import datetime
from bson.son import SON
from bson.objectid import ObjectId
from biocluster.api.database.base import Base, report_check


class BsaQc(Base):
    def __init__(self):
        super(BsaQc, self).__init__()
        """
        laste modified by hd 修改了mongo的连接20180222
        """
        self._project_type = "bsa"

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
        coll = self.db["sg_task"]
        main_id = coll.insert_one(insert_data).inserted_id
        print "add_sg_task success!"
        return main_id

    def add_sg_specimen(self, project_sn, task_id, old_name, new_name, clean_path, desc=None):
        for f in clean_path.split(";"):
            if not os.path.exists(f):
                raise Exception("{}所指定的路径不存在，请检查!".format(f))
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "old_name": old_name,
            "new_name": new_name,
            "clean_path": clean_path,
            "desc": desc if desc else "sg_specimen"
        }
        coll = self.db["sg_specimen"]
        main_id = coll.insert_one(insert_data).inserted_id
        print "add_sg_specimen success!"
        return main_id

    def add_sg_specimen_qc(self, task_id, specimen_id, qc_stat):
        """
        """
        if not os.path.exists(qc_stat):
            raise Exception("{}所指定的路径不存在，请检查".format(qc_stat))
        with open(qc_stat, "r") as f:
            for line in f:
                if line.startswith("#"):
                    pass
                else:
                    item = line.strip().split("\t")
                    insert_data = {
                        "task_id": task_id,
                        "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        "specimen_id": str(specimen_id),
                        "clean_reads": int(item[2]),
                        "clean_base": int(item[3]),
                        "gc_rate": float(item[9]),
                        "q30_rate": float(item[10])
                    }
                    break
        coll = self.db["sg_specimen_qc"]
        main_id = coll.insert_one(insert_data).inserted_id
        coll.update_one({"_id": main_id}, {"$set": {"main_id": main_id}})
        print "add_sg_specimen_qc success!"
        return main_id

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

    def add_sg_base_content(self, task_id, specimen_qc_id, qual_path):
        if not isinstance(specimen_qc_id, ObjectId):
            if isinstance(specimen_qc_id, types.StringTypes):
                specimen_qc_id = ObjectId(specimen_qc_id)
            else:
                raise Exception("specimen_qc_id必须为ObjectID或其对应的字符串")
        if not os.path.exists(qual_path):
            raise Exception("{}文件路径不存在，请检查!".format(qual_path))
        categories = []
        a_list, t_list, g_list, c_list, n_list = [], [], [], [], []
        sample = os.path.basename(qual_path).split(".")[0]
        with open(qual_path, "r") as f:
            header = f.readline()
            for line in f:
                item = line.strip().split("\t")
                categories.append(item[0])
                a_list.append(float(item[1]))
                t_list.append(float(item[2]))
                g_list.append(float(item[3]))
                c_list.append(float(item[4]))
                n_list.append(float(item[5]))
        main_id = self.add_sg_curve(task_id=task_id, specimen_qc_id=specimen_qc_id, name="碱基含量分布图",\
                  categories=categories, location="Bases content along raw reads", sample=sample, types=1)
        self.add_sg_curve_detail(curve_id=main_id, name="A", value=a_list)
        self.add_sg_curve_detail(curve_id=main_id, name="T", value=t_list)
        self.add_sg_curve_detail(curve_id=main_id, name="C", value=c_list)
        self.add_sg_curve_detail(curve_id=main_id, name="G", value=g_list)
        self.add_sg_curve_detail(curve_id=main_id, name="N", value=n_list)

    def add_sg_base_error(self, task_id, specimen_qc_id, qual_path):
        if not isinstance(specimen_qc_id, ObjectId):
            if isinstance(specimen_qc_id, types.StringTypes):
                specimen_qc_id = ObjectId(specimen_qc_id)
            else:
                raise Exception("specimen_qc_id必须为ObjectID或其对应的字符串")
        if not os.path.exists(qual_path):
            raise Exception("{}文件路径不存在，请检查!".format(qual_path))
        categories = []
        e_list = []
        sample = os.path.basename(qual_path).split(".")[0]
        with open(qual_path, "r") as f:
            header = f.readline()
            for line in f:
                item = line.strip().split("\t")
                categories.append(item[0])
                error = 1 / (10 ** (float(item[-1]) / 10))
                e_list.append(error)
        main_id = self.add_sg_curve(task_id=task_id, specimen_qc_id=specimen_qc_id, name="碱基错误率分布图",\
                  categories=categories, location="Mean error distribution along raw reads", sample=sample, types=1)
        self.add_sg_curve_detail(curve_id=main_id, name="", value=e_list)


# if __name__ == "__main__":
#     project_sn = "bsa_test"
#     task_id = "bsa_test"
#     member_id = "m_test"
#     member_type = 0
#     cmd_id = 0
#     a = BsaQc()
#     a.add_sg_task(member_id=member_id, member_type=member_type, cmd_id=cmd_id, project_sn=project_sn, task_id=task_id)
#     clean_path = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/01.fastq_qc/stat/B23XC-1.stat;/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/01.fastq_qc/stat/B23XC-1.stat"
#     qc_stat = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/01.fastq_qc/stat/B23XC-1.stat"
#     atgc_path = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/01.fastq_qc/atgc/B23XC-1.atgc"
#     qual_path = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/01.fastq_qc/qual/B23XC-1.qual"
#     specimen_id = a.add_sg_specimen(project_sn=project_sn, task_id=task_id, old_name="B23XC", new_name="B23XC-1", clean_path=clean_path, desc=None)
#     specimen_qc_id = a.add_sg_specimen_qc(task_id=task_id, specimen_id=specimen_id, qc_stat=qc_stat)
#     a.add_sg_base_content(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=atgc_path)
#     a.add_sg_base_error(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=qual_path)
#     qc_stat = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/01.fastq_qc/stat/C10XC-1.stat"
#     atgc_path = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/01.fastq_qc/atgc/C10XC-1.atgc"
#     qual_path = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/01.fastq_qc/qual/C10XC-1.qual"
#     specimen_id = a.add_sg_specimen(project_sn=project_sn, task_id=task_id, old_name="C10XC", new_name="C10XC-1", clean_path=clean_path, desc=None)
#     specimen_qc_id = a.add_sg_specimen_qc(task_id=task_id, specimen_id=specimen_id, qc_stat=qc_stat)
#     a.add_sg_base_content(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=qual_path)
#     a.add_sg_base_error(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=atgc_path)
#     qc_stat = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/01.fastq_qc/stat/ZJU_co-1.stat"
#     atgc_path = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/01.fastq_qc/atgc/ZJU_co-1.atgc"
#     qual_path = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/01.fastq_qc/qual/ZJU_co-1.qual"
#     specimen_id = a.add_sg_specimen(project_sn=project_sn, task_id=task_id, old_name="ZJU_co", new_name="ZJU_co-1", clean_path=clean_path, desc=None)
#     specimen_qc_id = a.add_sg_specimen_qc(task_id=task_id, specimen_id=specimen_id, qc_stat=qc_stat)
#     a.add_sg_base_content(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=atgc_path)
#     a.add_sg_base_error(task_id=task_id, specimen_qc_id=specimen_qc_id, qual_path=qual_path)
