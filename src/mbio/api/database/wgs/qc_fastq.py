# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.05.02

from api_base import ApiBase
import datetime
import os
import re


class QcFastq(ApiBase):
    def __init__(self, bind_object):
        """
        用于module:qc_stat.py中，根据task_id查询mongo表sg_task/sg_specimen_other,得到样本信息和对应的fastq
        """
        super(QcFastq, self).__init__(bind_object)
        self._project_type = "dna_wgs"

    def get_specimen_fastq(self, task_id, fastq_dir=None):
        """
        """
        specimen_info, rawdata_info = {}, {}
        if not fastq_dir:
            result = self.col_find_one("sg_task", {"task_id": task_id})
            if not result:
                raise Exception("没有在sg_task表中找到task_id:{}".format(task_id))
            fastq_dir = result["raw_fastq_path"]
        self.check_exists(fastq_dir)
        results = self.col_find("sg_specimen_other", {"task_id": task_id, "selected": "1"})
        if results.count() == 0:  # 修改bug查找不到要报错 modified by hd 20180503, 当结果为空的时候，可以去看
            # 下表中selected是否是字符串
            raise Exception("没有在sg_specimen_other表中找到task_id:{}对应的样本".format(task_id))
        for result in results:
            initial_name = result["initial_name"]
            analysis_name = initial_name.split("-")[0]
            files = result["file_name"].split(",")
            if analysis_name not in specimen_info.keys():
                specimen_info[analysis_name] = []
            rawdata_info[initial_name] = {}
            rawdata_info[initial_name]["lib"] = result["library"]
            # for f in files:
            #     file = os.path.join(fastq_dir, f)
            #     if not os.path.exists(file):
            #         raise Exception("fastq文件:{}不存在，请检查".format(file))
            #     if re.search(r".*R1.*", f):
            #         rawdata_info[initial_name]["l"] = file
            #     if re.search(r".*R2.*", f):
            #         rawdata_info[initial_name]["r"] = file
            # if "l" not in rawdata_info[initial_name].keys():
            #     rawdata_info[initial_name]["l"] = os.path.join(fastq_dir, files[0])
            # if "r" not in rawdata_info[initial_name].keys():
            #     rawdata_info[initial_name]["r"] = os.path.join(fastq_dir, files[1])
            if files[0] == files[1]:
                raise Exception("样本:{}的左右两端fastq文件相同，请检查!".format(initial_name))
            fastq_l = os.path.join(fastq_dir, files[0])
            fastq_r = os.path.join(fastq_dir, files[1])
            self.check_exists(fastq_l)
            self.check_exists(fastq_r)
            rawdata_info[initial_name]["l"] = fastq_l
            rawdata_info[initial_name]["r"] = fastq_r
            specimen_info[analysis_name].append(initial_name)
        return specimen_info, rawdata_info


if __name__ == "__main__":
    a = QcFastq(None)
    task_id = "tsg_29684"
    fastq_dir = "/mnt/ilustre/tsanger-data/rerewrweset/files/m_188/188_5af3dc68d1df4/glycine_max_raw"
    a.get_specimen_fastq(task_id, fastq_dir)
