# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 2018.0706

from api_base import ApiBase
import datetime
import time
import os
import re
import math


class GmapBase(ApiBase):
    def __init__(self, bind_object):
        """
        gamp项目基础数据导表， 主要是导一些wgs的结果数据，如果wgs导表函数直接能用的话，就直接用wgs的导表
        """
        super(GmapBase, self).__init__(bind_object)
        self._project_type = "dna_gmap"
        self.project_sn = self.bind_object.sheet.project_sn
        self.task_id = self.bind_object.sheet.id

    def add_sg_task(self, member_id, member_type, cmd_id):
        """
        sg_task
        """
        data_list = []
        insert_data = {
            "project_sn": self.project_sn,
            "task_id": self.task_id,
            "member_id": member_id,
            "member_type": member_type,
            "cmd_id": cmd_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "project_type": "gmap"  # gmap项目导表测试
        }
        data_list.append(insert_data)
        self.col_insert_data("sg_task", data_list)
        print "sg_task导入任务{}成功".format(self.task_id)

    def add_sg_specimen(self, results, pid=None, mid=None):
        """
        sg_specimen
        pid:父本, mid:母本
        """
        infos = self.col_find("sg_specimen", {"task_id": self.task_id})
        if infos.count() != 0:
            self.db['sg_specimen'].remove({"task_id": self.task_id})
        data_list = []
        specimen_dict = {}
        specimen_ids = []
        for result in results:
            analysis_name = result["initial_name"].split("-")[0]
            if analysis_name not in specimen_dict.keys():
                specimen_ids.append(analysis_name)
                specimen_dict[analysis_name] = {"init": [], "raw": [], "lib": "", "run": []}
            specimen_dict[analysis_name]["init"].append(result["initial_name"])
            specimen_dict[analysis_name]["raw"].append(result["file_name"])
            specimen_dict[analysis_name]["lib"] = result["library"]
            specimen_dict[analysis_name]["run"].append(result["batch"])
        try:
            specimen_ids.sort(key=lambda i: int(re.findall("\d+", i)[0]))
        except:
            specimen_ids.sort()
        if mid in specimen_ids:
            specimen_ids.remove(mid)
            specimen_ids.insert(0, mid)
        if pid in specimen_ids:
            specimen_ids.remove(pid)
            specimen_ids.insert(0, pid)
        for s in specimen_ids:
            insert_data = {
                "project_sn": self.project_sn,
                "task_id": self.task_id,
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "initial_name": "|".join(specimen_dict[s]["init"]),
                "old_name": s,
                "new_name": s,
                "run": "|".join(specimen_dict[s]["run"]),
                "library": specimen_dict[s]["lib"],
                "raw_path": "|".join(specimen_dict[s]["raw"]),
                "clean_path": "no",
                "desc": ""
            }
            data_list.append(insert_data)
        self.col_insert_data("sg_specimen", data_list)
        print "sg_specimen导表成功"

    def add_sg_specimen_other(self, file_path, task_id):
        """
        fastq.txt文件
        Sample Initial Name	Sample Analysis Name	Batch	Library	File Name
        vv_96-2	vv_96	batch01	WGS	vv_96-2.lane5.clean_1.fastq.gz,vv_96-2.lane5.clean_2.fastq.gz
        vv_96-1	vv_96	batch01	WGS	vv_96-1.0810.clean_1.fastq.gz,vv_96-1.0810.clean_2.fastq.gz
        vv_38	vv_38	batch01	WGS	vv_38.0810.clean_1.fastq.gz,vv_38.0810.clean_2.fastq.gz
        :param file_path:
        :param task_id:
        :return:
        """
        data_list = []
        with open(file_path, "r") as r:
            data = r.readlines()[1:]
            for line in data:
                temp = line.strip().split("\t")
                insert_data = {
                    "task_id": task_id,
                    "initial_name": temp[0],
                    "analysis_name": temp[1],
                    "batch": temp[2],
                    "library": temp[3],
                    "file_name": temp[4],
                    "selected": "1"
                }
                try:
                    insert_data["other_task_id"] = task_id.split("_")[1]
                except:
                    insert_data["other_task_id"] = task_id
                data_list.append(insert_data)
        if data_list:
            self.col_insert_data("sg_specimen_other", data_list)

    def get_chrs(self, ref_chrlist):
        """
        获取染色体列表
        :param ref_chrlist:
        :return:
        """
        chr_list = []
        with open(ref_chrlist, 'r') as f:
            lines = f.readlines()
            for i in lines:
                chrs = i.strip().split('\t')
                if chrs[0] not in chr_list:
                    chr_list.append(chrs[0])
        return chr_list

    def get_child_id(self, task_id):
        """
        查找sg_specimen_child中的所有的样本，然后整理成sg_child_list格式导入，返回sg_child_list id
        :param task_id:
        :return:
        """
        sample_list = []
        data_list = []
        result = self.col_find("sg_specimen_child", {"task_id": task_id, "selected": 1})
        if result.count() != 0:
            for m in result:
                if m['analysis_name'] not in sample_list:
                    sample_list.append(m['analysis_name'])
            data_list.append({
                "task_id": task_id,
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "name": "child_{}".format(time.time()),
                "num": len(sample_list),
                "spcimen_ids": sample_list
            })
            if data_list:
                main_id = self.col_insert_data("sg_child_list", data_list)
                self.update_db_record("sg_child_list", {"_id": main_id}, {"main_id": main_id})
                return True, ",".join(sample_list), len(sample_list)
        else:
            return False, "", 0

if __name__ == "__main__":
    a = GmapBase(None)
    member_id = ""
    member_type = 1
    cmd_id = 1
    a.project_sn = 'bug_test'
    a.task_id = 'bug_test'
    a.add_sg_specimen("tsanger_31601", "P1", "P2")
